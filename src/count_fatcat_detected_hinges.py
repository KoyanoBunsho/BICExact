import os
import pandas as pd
import re
from tqdm import tqdm


def main():
    output_dir = "fatcat_simulation_result"
    simulation_file = "output_simulation_file_for_fatcat.csv"
    df = pd.read_csv(simulation_file)
    hinge_count_list = []
    hinge_info_list = []
    hinge_positions_list = []
    for _, row in tqdm(df.iterrows(), total=len(df)):
        res_file = f"{output_dir}/{os.path.basename(row['p_pdb'])}_{os.path.basename(row['q_pdb'])}.aln"
        hinge_positions = extract_hinge_positions(res_file)
        hinge_positions_list.append(hinge_positions)
        hinge_count = extract_hinge_number(res_file)
        hinge_count_list.append(hinge_count)
        hinge_info_list.append(row["hinge_file_path"])
    df["hinge_index"] = hinge_positions_list
    df["detected_hinge_count"] = hinge_count_list
    df["hinge_file_path"] = hinge_info_list
    df.to_csv("fatcat_simulation_hinge_count_result.csv", index=False)


def extract_hinge_number(file_path):
    pattern = re.compile(r"Twists\s+(\d+)")
    try:
        with open(file_path, "r") as file:
            content = file.read()
        match = pattern.search(content)
        if match:
            return int(match.group(1))
    except Exception as e:
        print(f"Failed to open or read {file_path}: {e}")
    return None


def extract_hinge_positions(file_path):
    twist_pattern = re.compile(r"Twists\s+(\d+)")
    valid_line_pattern = re.compile(r"^[\d\s]+$")
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()
        for line in lines:
            twist_match = twist_pattern.search(line)
            if twist_match:
                twists = int(twist_match.group(1))
                if twists > 0:
                    break
        else:
            return ""
        final_sequence = ""
        for line in lines:
            if valid_line_pattern.match(line.strip()):
                final_sequence += line.replace(" ", "").strip()
        return " : ".join(map(str, find_switch_indices(final_sequence)))
    except Exception as e:
        print(f"Failed to process {file_path}: {e}")
    return ""


def find_switch_indices(sequence):
    indices = []
    for i in range(1, len(sequence)):
        if sequence[i] != sequence[i - 1]:
            indices.append(i)
    return indices


if __name__ == "__main__":
    main()
