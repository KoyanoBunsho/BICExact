import os
import pandas as pd
from tqdm import tqdm
import re
import numpy as np


def main():
    simulation_file = (
        "output_simulation_file_for_fatcat.csv"  # Use the same file as fatacat
    )
    df = pd.read_csv(simulation_file)
    hinge_count_list = []
    hinge_info_list = []
    hinge_positions_list = []
    for _, row in tqdm(df.iterrows(), total=len(df)):
        pdb1_basename = os.path.basename(row["p_pdb"]).replace(".pdb", "")
        pdb2_basename = os.path.basename(row["q_pdb"]).replace(".pdb", "")
        res_file = f"{pdb1_basename}_{pdb2_basename}_info"
        hinge_positions = extract_hinge_positions(res_file)
        hinge_positions_list.append(hinge_positions)
        hinge_count = count_bending_residues(res_file)
        hinge_count_list.append(hinge_count)
        hinge_info_list.append(row["hinge_file_path"])
    df["hinge_index"] = hinge_positions_list
    df["detected_hinge_count"] = hinge_count_list
    df["hinge_file_path"] = hinge_info_list
    df.to_csv("dyndom_simulation_hinge_count_result.csv", index=False)


def count_bending_residues(file_path):
    pattern = re.compile(r"BENDING RESIDUES:\s*\d+\s*-\s*\d+")
    hinge_count = 0
    with open(file_path, "r") as file:
        for line in file:
            if pattern.search(line):
                hinge_count += 1
    return hinge_count


def extract_hinge_positions(file_path):
    pattern = re.compile(r"BENDING RESIDUES:\s*(\d+)\s*-\s*(\d+)")
    hinge_positions = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                match = pattern.search(line)
                if match:
                    start = int(match.group(1))
                    end = int(match.group(2))
                    center = int(np.ceil((end - start) / 2)) + start
                    hinge_positions.append(center)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    return " : ".join(map(str, hinge_positions))


if __name__ == "__main__":
    main()
