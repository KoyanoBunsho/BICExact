#ifndef BICEXACT_STRUCTURE_IO_H
#define BICEXACT_STRUCTURE_IO_H

#include "rmsd_struct.h"
#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace structure_io {
using Residue = std::pair<int, char>;
using Coordinates = std::map<Residue, Eigen::Vector3d>;
const std::size_t maxResidues = 10000;

inline std::string lower(std::string s) {
  for (auto &c : s)
    c = std::tolower(static_cast<unsigned char>(c));
  return s;
}

// Bound memory even for an input with no newlines. Both fixed-column PDB and
// CIF 1.1 records fit well below this limit.
class Lines {
  std::ifstream file;
  static constexpr std::size_t capacity = 1024 * 1024;
  std::unique_ptr<char[]> buffer{new char[capacity]};

public:
  explicit Lines(const std::string &path) : file(path) {
    if (!file)
      throw std::invalid_argument("Unable to open structure");
  }

  bool next(std::string &line) {
    if (!file.getline(buffer.get(), capacity)) {
      if (file.eof() && !file.bad() && !file.gcount())
        return false;
      throw std::invalid_argument("Unreadable or overlong structure record");
    }
    line.assign(buffer.get(), file.gcount() - (file.eof() ? 0 : 1));
    if (line.find('\0') != std::string::npos)
      throw std::invalid_argument("Invalid structure record");
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    return true;
  }
};

inline int residueNumber(const std::string &s) {
  std::size_t end;
  int n = std::stoi(s, &end);
  if (s.find_first_not_of(' ', end) != std::string::npos)
    throw std::invalid_argument("Invalid residue number");
  return n;
}

inline double coordinate(const std::string &s) {
  std::size_t end;
  double v = std::stod(s, &end);
  if (s.find_first_not_of(' ', end) != std::string::npos || !std::isfinite(v) ||
      std::abs(v) > 1e6)
    throw std::invalid_argument("Invalid atom coordinate");
  return v;
}

inline void add(Coordinates &coords, const std::string &seq, char insertion,
                const std::string &x, const std::string &y,
                const std::string &z) {
  if (insertion != ' ' && !std::isalnum(static_cast<unsigned char>(insertion)))
    throw std::invalid_argument("Invalid insertion code");
  coords.emplace(Residue{residueNumber(seq), insertion},
                 Eigen::Vector3d{coordinate(x), coordinate(y), coordinate(z)});
  if (coords.size() > maxResidues)
    throw std::invalid_argument("Too many residues");
}

inline void pdb(Lines &lines, const std::string &chain, Coordinates &coords) {
  std::string line;
  while (lines.next(line)) {
    if (line.compare(0, 6, "ENDMDL") == 0)
      break;
    if (line.compare(0, 4, "ATOM") != 0 || line.size() < 5 ||
        !std::isspace(static_cast<unsigned char>(line[4])))
      continue;
    if (line.compare(0, 6, "ATOM  ") == 0 && line.size() >= 54 &&
        line[11] == ' ' && line[20] == ' ' && line.compare(27, 3, "   ") == 0) {
      const auto atom = line.substr(12, 4);
      const auto start = atom.find_first_not_of(' ');
      if (chain.size() == 1 && chain[0] == line[21] &&
          start != std::string::npos &&
          atom.compare(start, atom.find_last_not_of(' ') - start + 1, "CA") ==
              0)
        add(coords, line.substr(22, 4), line[26], line.substr(30, 8),
            line.substr(38, 8), line.substr(46, 8));
    } else {
      // Earlier releases also accepted whitespace-separated ATOM records.
      std::istringstream fields(line);
      std::string ignored, atom, asym, seq, x, y, z;
      if (!(fields >> ignored >> ignored >> atom >> ignored >> asym >> seq >>
            x >> y >> z))
        throw std::invalid_argument("Incomplete PDB record");
      if (atom == "CA" && asym == chain)
        add(coords, seq, ' ', x, y, z);
    }
  }
}

struct Token {
  std::string value;
  bool quoted = false;
};

class CifTokens {
  Lines &lines;
  std::string line;
  std::size_t pos = 0;

public:
  explicit CifTokens(Lines &source) : lines(source) {}

  bool next(Token &token) {
    for (;;) {
      pos = line.find_first_not_of(" \t\r", pos);
      if (pos != std::string::npos && line[pos] != '#')
        break;
      if (!lines.next(line))
        return false;
      pos = 0;
    }
    const char first = line[pos];
    token.quoted = true;
    if (pos == 0 && first == ';') {
      // Text fields occur in metadata. They cannot be numeric atom fields.
      while (lines.next(line)) {
        if (!line.empty() && line.front() == ';') {
          pos = 1;
          token.value.clear();
          return true;
        }
      }
      throw std::invalid_argument("Unterminated CIF text");
    }
    if (first == '"' || first == '\'') {
      const auto start = ++pos;
      while ((pos = line.find(first, pos)) != std::string::npos) {
        ++pos;
        if (pos == line.size() ||
            std::isspace(static_cast<unsigned char>(line[pos]))) {
          token.value = line.substr(start, pos - start - 1);
          return true;
        }
      }
      throw std::invalid_argument("Unterminated CIF quote");
    }
    token.quoted = false;
    auto start = pos;
    pos = line.find_first_of(" \t\r", start);
    token.value = line.substr(start, pos - start);
    return true;
  }

  void legacyRow(const std::string &chain, Coordinates &coords) {
    std::istringstream fields(line.substr(pos));
    std::string ignored, atom, asym, seq, insertion, x, y, z;
    if (!(fields >> ignored >> ignored >> atom >> ignored >> ignored >> asym >>
          ignored >> seq >> insertion >> x >> y >> z))
      throw std::invalid_argument("Incomplete CIF atom row");
    if (atom == "CA" && asym == chain)
      add(coords, seq, ' ', x, y, z);
    pos = line.size();
  }
};

inline bool cifControl(const Token &t) {
  if (t.quoted || t.value.empty())
    return false;
  if (t.value[0] == '_')
    return true;
  if (t.value.find('_') == std::string::npos)
    return false;
  const auto s = lower(t.value);
  return s == "loop_" || s == "stop_" || s == "global_" ||
         s.compare(0, 5, "data_") == 0 || s.compare(0, 5, "save_") == 0;
}

struct AtomSite {
  std::string group, atom, chain, seq, authAtom, authChain, authSeq, insertion;
  std::string x, y, z, model;
  using Member = std::string AtomSite::*;

  static Member member(const std::string &tag) {
    static const std::map<std::string, Member> fields{
        {"_atom_site.group_pdb", &AtomSite::group},
        {"_atom_site.label_atom_id", &AtomSite::atom},
        {"_atom_site.label_asym_id", &AtomSite::chain},
        {"_atom_site.label_seq_id", &AtomSite::seq},
        {"_atom_site.auth_atom_id", &AtomSite::authAtom},
        {"_atom_site.auth_asym_id", &AtomSite::authChain},
        {"_atom_site.auth_seq_id", &AtomSite::authSeq},
        {"_atom_site.pdbx_pdb_ins_code", &AtomSite::insertion},
        {"_atom_site.cartn_x", &AtomSite::x},
        {"_atom_site.cartn_y", &AtomSite::y},
        {"_atom_site.cartn_z", &AtomSite::z},
        {"_atom_site.pdbx_pdb_model_num", &AtomSite::model}};
    auto it = fields.find(tag);
    return it == fields.end() ? nullptr : it->second;
  }

  static bool present(const std::string &v) {
    return !v.empty() && v != "." && v != "?";
  }

  void collect(const std::string &selectedChain, bool pdbNumbering,
               std::string &firstModel, Coordinates &coords) const {
    if (firstModel.empty())
      firstModel = model;
    if (model != firstModel || group != "ATOM" ||
        (present(atom) ? atom : authAtom) != "CA")
      return;
    // Mixed-format comparisons use author identifiers to match the PDB file.
    const bool useAuthChain =
        present(authChain) && (pdbNumbering || !present(chain));
    if ((useAuthChain ? authChain : chain) != selectedChain)
      return;
    const bool useAuthSeq = present(authSeq) && (pdbNumbering || !present(seq));
    char code = ' ';
    if (useAuthSeq && present(insertion)) {
      if (insertion.size() != 1)
        throw std::invalid_argument("Invalid insertion code");
      code = insertion[0];
    }
    add(coords, useAuthSeq ? authSeq : seq, code, x, y, z);
  }
};

inline void cif(Lines &lines, const std::string &chain, bool pdbNumbering,
                Coordinates &coords) {
  CifTokens tokens(lines);
  Token t;
  bool more = tokens.next(t);
  if (more && !t.quoted && (t.value == "ATOM" || t.value == "HETATM")) {
    do {
      if (!t.quoted && t.value == "ATOM")
        tokens.legacyRow(chain, coords);
    } while (tokens.next(t));
    return;
  }
  AtomSite single;
  std::string firstModel;
  while (more) {
    const auto tag = t.quoted ? "" : lower(t.value);
    if (tag == "loop_") {
      AtomSite row;
      std::size_t width = 0;
      // Keep only recognized atom fields, independently of the number of tags.
      std::map<std::size_t, AtomSite::Member> columns;
      while ((more = tokens.next(t)) && !t.quoted && !t.value.empty() &&
             t.value.front() == '_') {
        if (auto member = AtomSite::member(lower(t.value))) {
          if (std::any_of(columns.begin(), columns.end(),
                          [member](const auto &col) { return col.second == member; }))
            throw std::invalid_argument("Duplicate CIF atom field");
          columns.emplace(width, member);
        }
        ++width;
      }
      if (!width)
        throw std::invalid_argument("Missing CIF loop columns");
      std::size_t i = 0;
      while (more && !cifControl(t)) {
        auto col = columns.find(i);
        if (col != columns.end())
          (row.*(col->second)).swap(t.value);
        if (++i == width) {
          if (!columns.empty())
            row.collect(chain, pdbNumbering, firstModel, coords);
          i = 0;
        }
        more = tokens.next(t);
      }
      if (i)
        throw std::invalid_argument("Incomplete CIF loop row");
      if (!columns.empty())
        return;
    } else if (!tag.empty() && tag.front() == '_') {
      if (!tokens.next(t) || cifControl(t))
        throw std::invalid_argument("Missing CIF value");
      if (auto member = AtomSite::member(tag))
        single.*member = t.value;
      more = tokens.next(t);
    } else {
      if (tag.compare(0, 5, "data_") == 0 && !single.group.empty())
        break;
      more = tokens.next(t);
    }
  }
  single.collect(chain, pdbNumbering, firstModel, coords);
}

inline std::string extension(const std::string &path) {
  auto dot = path.rfind('.');
  return dot == std::string::npos ? "" : lower(path.substr(dot));
}

inline Coordinates read(const std::string &path, const std::string &chain,
                        bool pdbNumbering) {
  Lines lines(path);
  Coordinates coords;
  auto ext = extension(path);
  if (ext == ".cif")
    cif(lines, chain, pdbNumbering, coords);
  else if (ext == ".pdb" || ext == ".ent")
    pdb(lines, chain, coords);
  else
    throw std::invalid_argument("Unsupported structure format");
  return coords;
}
} // namespace structure_io
#endif
