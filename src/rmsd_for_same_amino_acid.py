import os
import pandas as pd
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed
from tqdm import tqdm
import rmsd


def main():
    pdb_chain_pairs = pd.read_csv("pdb_chain_id_pair_for_rmsd.csv")
    pdb_directory = "all_pdb"
    n_jobs = os.cpu_count() or 1
    print(f"{n_jobs} parallel jobs running")

    rmsd_result_list = Parallel(n_jobs=n_jobs)(
        delayed(process_pair)(row, pdb_directory)
        for _, row in tqdm(pdb_chain_pairs.iterrows())
    )

    pd.DataFrame(
        [item for sublist in rmsd_result_list for item in sublist if item]
    ).to_csv("rmsd_distribution_240502.csv", index=False)


def extract_ca_coordinates(pdb_id, chain_id, pdb_directory):
    file_path = os.path.join(pdb_directory, f"pdb{pdb_id}.ent.gz")
    ppdb = PandasPdb().read_pdb(file_path)
    df = ppdb.df["ATOM"]
    ca_df = df[
        (df["chain_id"] == chain_id) & (df["atom_name"] == "CA")
    ].drop_duplicates(subset=["residue_number"])
    return ca_df[["x_coord", "y_coord", "z_coord"]].values


def calculate_rmsd(coords1, coords2):
    if coords1.shape[0] == coords2.shape[0]:
        coords1_centered = coords1 - rmsd.centroid(coords1)
        coords2_centered = coords2 - rmsd.centroid(coords2)
        return rmsd.kabsch_rmsd(coords1_centered, coords2_centered)
    else:
        return None


def process_pair(row, pdb_directory):
    p_pdb, q_pdb = row["p_pdb"], row["q_pdb"]
    p_pdb_id, p_chain_id = p_pdb.split("_")
    q_pdb_id, q_chain_id = q_pdb.split("_")

    p_coordinates = extract_ca_coordinates(p_pdb_id, p_chain_id, pdb_directory)
    q_coordinates = extract_ca_coordinates(q_pdb_id, q_chain_id, pdb_directory)

    if p_coordinates.size and q_coordinates.size:
        rmsd_result = calculate_rmsd(p_coordinates, q_coordinates)
        if rmsd_result is not None:
            return [{"p_pdb": p_pdb, "q_pdb": q_pdb, "RMSD": rmsd_result}]
    return [None]


if __name__ == "__main__":
    main()
