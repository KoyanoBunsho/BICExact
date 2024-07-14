import math
import os
from copy import deepcopy
from glob import glob

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed
from scipy.spatial.transform import Rotation as R
from scipy.stats import truncnorm
from tqdm import tqdm

p = pd.read_csv("b_factors_ratios_strict.csv")["Flexible"]


def main():
    # pdb_files = glob("all_pdb/*.ent.gz")
    pdb_files = pd.read_csv("selected_files.csv", header=None)[0].to_list()
    os.makedirs("simulation_data", exist_ok=True)
    os.makedirs("simulation_data_info", exist_ok=True)
    n_jobs = os.cpu_count()
    print(f"{n_jobs} parallel")
    simulation_data_list = Parallel(n_jobs=n_jobs)(
        delayed(process_pdb_file)(pdb_file)
        for pdb_file in tqdm(pdb_files, total=len(pdb_files))
    )
    simulation_data_flat_list = [
        item for sublist in simulation_data_list for item in sublist
    ]
    sigma = 1.55 / np.sqrt(3)
    pd.DataFrame(simulation_data_flat_list).to_csv(
        f"simulation_data_sigma_{sigma}.csv", index=False
    )


def process_pdb_file(pdb_file):
    ppdb, df, chain_id = load_and_filter_pdb(pdb_file)
    if ppdb is None or df is None or chain_id is None:
        print(f"Skipping {pdb_file} due to missing CA atoms.")
        return []
    pdb_id = pdb_file.split("/")[1].split(".")[0]
    threshold = calc_threshold(df.copy())
    sigma = 1.55 / np.sqrt(3)
    np.random.seed()
    simulation_data_list = []
    for hinge_num in range(min(6, df["unique_residue_number"].max() // 21)):
        rot_df = df.copy()
        cnt = 0
        while True:
            if cnt > 10:
                break
            cnt += 1
            if hinge_num > 0:
                hinge_indices = select_hinge_indices(
                    2, df["unique_residue_number"].max() - 20, hinge_num, selected=[]
                )
                if len(hinge_indices) < hinge_num:
                    continue
                hinge_indices.sort()
                if not_valid_hinge_indices(hinge_indices=hinge_indices):
                    continue
                for hinge_index in hinge_indices:
                    angle1 = np.random.uniform(low=0, high=360, size=1)[0]
                    angle2 = np.random.uniform(low=0, high=360, size=1)[0]
                    u = np.random.uniform(low=0, high=1, size=1)[0]
                    v = np.random.uniform(low=0, high=1, size=1)[0]
                    rot_df = random_rotate_residues(rot_df, hinge_index, u, v)
                    if len(rot_df) == 0:
                        break
                rot_df = add_noise_to_df(rot_df, sigma)

            else:
                angle1 = None
                angle2 = None
                hinge_indices = None
                rot_df = add_noise_to_df(rot_df, sigma)
            if save_structure_if_valid(
                deepcopy(ppdb), rot_df, pdb_id, threshold, hinge_num, chain_id, sigma
            ):
                if hinge_num > 0:
                    simulation_data_list += [
                        {
                            "pdb_id": pdb_id,
                            "chain_id": chain_id,
                            "angle1": angle1,
                            "angle2": angle2,
                            "hinge_indices": " : ".join(
                                [str(hinge_index) for hinge_index in hinge_indices]
                            ),
                            "hinge_num": hinge_num,
                        }
                    ]
                    pd.DataFrame(
                        {
                            "pdb_id": pdb_id,
                            "chain_id": chain_id,
                            "angle1": angle1,
                            "angle2": angle2,
                            "hinge_indices": " : ".join(
                                [str(hinge_index) for hinge_index in hinge_indices]
                            ),
                            "hinge_num": hinge_num,
                        },
                        index=[0],
                    ).to_csv(
                        f"simulation_data_info/{pdb_id}_{chain_id}_hinge_{hinge_num}_sigma{sigma}.csv",
                        index=False,
                    )
                else:
                    simulation_data_list += [
                        {
                            "pdb_id": pdb_id,
                            "chain_id": chain_id,
                            "angle1": angle1,
                            "angle2": angle2,
                            "hinge_indices": hinge_indices,
                            "hinge_num": hinge_num,
                        }
                    ]
                break
            else:
                rot_df = df.copy()
    return simulation_data_list


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


def not_valid_hinge_indices(hinge_indices):
    if len(hinge_indices) > 1:
        for i in range(len(hinge_indices) - 1):
            if hinge_indices[i + 1] - hinge_indices[i] < 5:
                return True
        return False
    else:
        return False


def calc_threshold(df):
    threshold = np.inf
    df = df[df["atom_name"] == "CA"]
    coords = df[["x_coord", "y_coord", "z_coord"]].values
    for i in range(1, len(coords) - 1):
        for j in range(i):
            dist = np.sqrt(np.dot(coords[i] - coords[j], coords[i] - coords[j]))
            threshold = min(threshold, dist)
    return threshold


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
        ppdb_copy.to_pdb(f"simulation_data/{pdb_id}_{first_chain_id}_original.pdb")
        codes, _ = pd.factorize(filtered_df["residue_number"], sort=False)
        filtered_df["unique_residue_number"] = codes + 1
        return ppdb, filtered_df, first_chain_id
    else:
        return None, None, None


def random_rotate_residues(df, hinge_index, u, v):
    if len(df) == 0:
        print("case 1")
        return pd.DataFrame()
    if df[
        (df["unique_residue_number"] == hinge_index) & (df["atom_name"] == "CA")
    ].empty:
        print("case 2")
        return pd.DataFrame()
    filtered_df = df[
        (df["unique_residue_number"] == hinge_index) & (df["atom_name"] == "CA")
    ][["x_coord", "y_coord", "z_coord"]]
    if len(filtered_df) > 0:
        pivot = filtered_df.values[0]
    else:
        return pd.DataFrame()
    if df[
        (df["unique_residue_number"] == hinge_index + 1) & (df["atom_name"] == "CA")
    ].empty:
        return pd.DataFrame()
    old_point = df[
        (df["unique_residue_number"] == hinge_index + 1) & (df["atom_name"] == "CA")
    ][["x_coord", "y_coord", "z_coord"]].values[0]
    rot_radius = np.linalg.norm(old_point - pivot)
    z = -2 * u + 1
    x = rot_radius * (np.sqrt(1 - z * z) * np.cos(2 * np.pi * v))
    y = rot_radius * (np.sqrt(1 - z * z) * np.sin(2 * np.pi * v))
    z *= rot_radius
    new_point = np.array([x, y, z])
    rot_axis = np.cross(old_point - pivot, new_point)
    cos_theta = np.dot(new_point, old_point - pivot) / (
        np.linalg.norm(new_point) * np.linalg.norm(old_point - pivot)
    )
    theta = np.arccos(cos_theta)
    rotation_matrix = R.from_rotvec(
        rot_axis * theta / np.linalg.norm(rot_axis)
    ).as_matrix()
    try:
        assert np.isclose(
            np.dot(rotation_matrix, old_point - pivot) + pivot, new_point + pivot
        ).all()
    except AssertionError:
        print("Assertion failed.")
        print("Rotation Matrix:\n", rotation_matrix)
        print("Old point + pivot:", np.dot(rotation_matrix, old_point - pivot) + pivot)
        print("New point + pivot:", new_point + pivot)
        raise
    for index, row in df[df["unique_residue_number"] > hinge_index].iterrows():
        coord = np.array([row["x_coord"], row["y_coord"], row["z_coord"]])
        rotated_coord = np.dot(rotation_matrix, (coord - pivot).T) + pivot.T
        df.at[index, "x_coord"] = rotated_coord[0]
        df.at[index, "y_coord"] = rotated_coord[1]
        df.at[index, "z_coord"] = rotated_coord[2]
    return df


def test_not_collision_efficient(df, threshold):
    if len(df) == 0:
        return False
    P = []
    df = df[df["atom_name"] == "CA"]
    for coord in df[["x_coord", "y_coord", "z_coord"]].values:
        P.append(Point(x=coord[0], y=coord[1], z=coord[2]))
    n = len(P)
    min_dist = closest(P, n)
    # print(min_dist)
    return min_dist > threshold


def add_noise_to_df(rot_df, sigma):
    for index in rot_df.index:
        rot_df.at[index, "x_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
        rot_df.at[index, "y_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
        rot_df.at[index, "z_coord"] += np.random.normal(loc=0.0, scale=sigma, size=1)[0]
    return rot_df


def add_norm_noise_to_df(rot_df, sigma):
    for index in rot_df.index:
        epsilon = generate_epsilon(sigma=sigma)
        rot_df.at[index, "x_coord"] += epsilon[0]
        rot_df.at[index, "y_coord"] += epsilon[1]
        rot_df.at[index, "z_coord"] += epsilon[2]
    return rot_df


def save_structure_if_valid(ppdb, df, pdb_id, threshold, hinge_num, chain_id, sigma):
    if test_not_collision_efficient(df.copy(), threshold=threshold):
        output_filename = (
            f"simulation_data/{pdb_id}_{chain_id}_hinge_{hinge_num}_sigma{sigma}.pdb"
        )
        ppdb.df["ATOM"] = df
        ppdb.to_pdb(output_filename)
        return True
    else:
        # print("Collision")
        return False


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def dist(p1, p2):
    return math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)


def bruteForce(P, n):
    min_dist = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            if dist(P[i], P[j]) < min_dist:
                min_dist = dist(P[i], P[j])
    return min_dist


def stripClosest(strip, size, d):
    min_dist = d
    strip = sorted(strip, key=lambda point: point.y)

    for i in range(size):
        for j in range(i + 1, size):
            if (strip[j].y - strip[i].y) >= min_dist:
                break
            if dist(strip[i], strip[j]) < min_dist:
                min_dist = dist(strip[i], strip[j])
    return min_dist


def closestUtil(P, n):
    if n <= 3:
        return bruteForce(P, n)
    mid = n // 2
    midPoint = P[mid]
    dl = closestUtil(P[:mid], mid)
    dr = closestUtil(P[mid:], n - mid)
    d = min(dl, dr)
    strip = []
    for i in range(n):
        if abs(P[i].x - midPoint.x) < d:
            strip.append(P[i])
    return min(d, stripClosest(strip, len(strip), d))


def closest(P, n):
    P = sorted(P, key=lambda point: point.x)
    return closestUtil(P, n)


def sample_truncated_normal(mean=0, sd=1, lower=0, upper=10, size=1):
    return truncnorm((lower - mean) / sd, (upper - mean) / sd, loc=mean, scale=sd).rvs(
        size
    )


def generate_epsilon(mean=0, sigma=1, lower=0, upper=np.inf):
    u = np.random.uniform(low=0, high=1, size=1)[0]
    v = np.random.uniform(low=0, high=1, size=1)[0]
    z = -2 * u + 1
    x = np.sqrt(1 - z * z) * np.cos(2 * np.pi * v)
    y = np.sqrt(1 - z * z) * np.sin(2 * np.pi * v)
    norm = sample_truncated_normal(mean, sigma, lower, upper)
    epsilon = np.array([x, y, z]) * norm
    return epsilon


if __name__ == "__main__":
    main()
