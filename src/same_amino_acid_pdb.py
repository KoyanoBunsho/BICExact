from collections import defaultdict
from glob import glob
import pandas as pd
from biopandas.pdb import PandasPdb
from tqdm import tqdm
from joblib import Parallel, delayed
import os

amino_acid_dict = {
    "ARG": "R",
    "HIS": "H",
    "LYS": "K",
    "ASP": "D",
    "GLU": "E",
    "SER": "S",
    "THR": "T",
    "ASN": "N",
    "GLN": "Q",
    "CYS": "C",
    "GLY": "G",
    "PRO": "P",
    "ALA": "A",
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "MET": "M",
    "PHE": "F",
    "TYR": "Y",
    "TRP": "W",
}


def main():
    similar_pdb_files = find_pdb_files_with_similar_sequences("all_pdb")
    pdb_chain_id_pair = []
    for _, files in similar_pdb_files.items():
        for i in range(1, len(files)):
            for j in range(i):
                pdb_chain_id_pair += [{"p_pdb": files[i], "q_pdb": files[j]}]
    pd.DataFrame(pdb_chain_id_pair).to_csv(
        "pdb_chain_id_pair_for_rmsd.csv", index=False
    )


def read_pdb_file(file_path):
    ppdb = PandasPdb().read_pdb(file_path)
    df = ppdb.df["ATOM"]
    dna_rna_residues = {"DA", "DC", "DG", "DT"}
    if df["residue_name"].isin(dna_rna_residues).any():
        return "", ""

    sequences = []
    chain_ids = df["chain_id"].unique()

    for chain_id in chain_ids:
        df_chain = df[
            (df["chain_id"] == chain_id) & (df["atom_name"] == "CA")
        ].drop_duplicates(subset=["residue_number"])
        amino_acid_sequence = "".join(
            df_chain["residue_name"].map(amino_acid_dict).fillna("")
        )
        if amino_acid_sequence:
            sequences.append((amino_acid_sequence, chain_id))
    return sequences


def process_file(file):
    sequences_chain = read_pdb_file(file)
    if sequences_chain != ("", ""):
        pdb_id = file.split("/")[-1].split(".")[0][3:]
        return [
            (amino_acid_sequence, f"{pdb_id}_{use_chain_id}")
            for amino_acid_sequence, use_chain_id in sequences_chain
        ]
    return []


def find_pdb_files_with_similar_sequences(directory):
    files = glob(f"{directory}/*.ent.gz")
    sequences = defaultdict(list)

    results = Parallel(n_jobs=os.cpu_count())(
        delayed(process_file)(file) for file in tqdm(files)
    )

    for result in results:
        for amino_acid_sequence, pdb_chain_id in result:
            sequences[amino_acid_sequence].append(pdb_chain_id)

    similar_sequences = {
        seq: paths for seq, paths in sequences.items() if len(paths) > 1
    }
    return similar_sequences


if __name__ == "__main__":
    main()
