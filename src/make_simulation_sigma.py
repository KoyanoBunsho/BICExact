import os
from copy import deepcopy

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
import rmsd
from joblib import Parallel, delayed
from tqdm import tqdm


def main():
    # pdb_files = glob("all_pdb/*.ent.gz")
    pdb_files = pd.read_csv("selected_files.csv", header=None)[0].to_list()
    n_jobs = os.cpu_count()
    print(f"{n_jobs} parallel")
    sigma_list = [0.5, 1.0, 1.5, 1.55, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    var_list = []
    for sigma in sigma_list:
        sigma /= np.sqrt(3)
        simulation_data_list = Parallel(n_jobs=n_jobs)(
            delayed(process_pdb_file)(pdb_file, sigma)
            for pdb_file in tqdm(pdb_files, total=len(pdb_files))
        )
        simulation_data_flat_list = [
            item for sublist in simulation_data_list for item in sublist
        ]
        res = 0
        for rmsd_res in simulation_data_flat_list:
            res += rmsd_res * rmsd_res
        res /= len(simulation_data_flat_list)
        res = np.sqrt(res)
        var_list.append(res)
    pd.DataFrame({"sigma": sigma_list, "variance": var_list}).to_csv(
        "sigma_variance_relationship.csv", index=False
    )


def process_pdb_file(pdb_file, sigma):
    ppdb, df, chain_id = load_and_filter_pdb(pdb_file)
    if ppdb is None or df is None or chain_id is None:
        print(f"Skipping {pdb_file} due to missing CA atoms.")
        return []
    np.random.seed()
    rot_df = add_noise_to_df(df.copy(), sigma)
    df_ca = df[df["atom_name"] == "CA"]
    rot_df_ca = rot_df[rot_df["atom_name"] == "CA"]
    coords1 = df_ca[["x_coord", "y_coord", "z_coord"]].to_numpy()
    coords2 = rot_df_ca[["x_coord", "y_coord", "z_coord"]].to_numpy()
    coords1_centered = coords1 - rmsd.centroid(coords1)
    coords2_centered = coords2 - rmsd.centroid(coords2)
    rmsd_result = rmsd.kabsch_rmsd(coords1_centered, coords2_centered)
    return [rmsd_result]


def select_hinge_indices(start, end, num_hinges, selected=[]):
    if num_hinges == 0:
        return selected
    elif len(selected) == 0:
        if start > end:
            return selected
        chosen_index = np.random.randint(start, end + 1)
        selected.append(chosen_index)
        return select_hinge_indices(start, end, num_hinges - 1, selected)
    else:
        excluded_range = set(
            range(max(start, selected[-1] - 20), min(end, selected[-1] + 20) + 1)
        )
        available = [
            i
            for i in range(start, end + 1)
            if i not in excluded_range and i not in selected
        ]
        if len(available) == 0 or len(available) < num_hinges - len(selected):
            return selected
        chosen_index = np.random.choice(available)
        selected.append(chosen_index)
        return select_hinge_indices(start, end, num_hinges - 1, selected)


def load_and_filter_pdb(filename):
    ppdb = PandasPdb().read_pdb(filename)
    if ppdb.df["ATOM"]["chain_id"].empty:
        return None, None, None
    first_chain_id = ppdb.df["ATOM"]["chain_id"].iloc[0]
    filtered_df = ppdb.df["ATOM"][ppdb.df["ATOM"]["chain_id"] == first_chain_id]
    if not filtered_df[filtered_df["atom_name"] == "CA"].empty:
        ppdb_copy = deepcopy(ppdb)
        ppdb_copy.df["ATOM"] = filtered_df
        pdb_id = filename.split("/")[-1].split(".")[0]
        codes, _ = pd.factorize(filtered_df["residue_number"], sort=False)
        filtered_df["unique_residue_number"] = codes + 1
        return ppdb, filtered_df, first_chain_id
    else:
        return None, None, None


def add_noise_to_df(rot_df, sigma):
    for index in rot_df.index:
        rot_df.at[index, "x_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
        rot_df.at[index, "y_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
        rot_df.at[index, "z_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
    return rot_df


if __name__ == "__main__":
    main()
