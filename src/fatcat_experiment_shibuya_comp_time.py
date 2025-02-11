import pandas as pd
import subprocess
import os
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed  # Joblibをインポート


def main():
    csv_file = "pdb_12_with_hinges_lower.csv"
    os.makedirs("dyndom_pdb_shibuya", exist_ok=True)
    if not os.path.exists(csv_file):
        print("Error: CSV file does not exist.")
        return
    pdb_pairs = load_pdb_pairs(csv_file)
    os.chdir("all_pdb")
    os.makedirs("fatcat_result_shibuya", exist_ok=True)

    num_cores = os.cpu_count()
    Parallel(n_jobs=num_cores)(
        delayed(extract_pdb)(pair) for index, pair in pdb_pairs.iterrows()
    )
    results = Parallel(n_jobs=num_cores)(
        delayed(run_fatcat)(pair) for index, pair in pdb_pairs.iterrows()
    )
    execution_time_list = [
        {"p_pdb": pair["p_pdb"], "q_pdb": pair["q_pdb"], "execution_time": duration}
        for pair, duration in zip(pdb_pairs.to_dict("records"), results)
    ]
    pd.DataFrame(execution_time_list).to_csv(
        "fatcat_execution_time_shibuya_improved.csv", index=False
    )


def load_pdb_pairs(csv_file):
    data = pd.read_csv(csv_file)
    return data


def run_fatcat(pair):
    pdb1 = pair["p_pdb"]
    pdb2 = pair["q_pdb"]
    output_filename = f"fatcat_result_shibuya/{pdb1}_{pdb2}"
    command = [
        "./FATCAT_speed",
        "-p1",
        f"{pdb1}.pdb",
        "-p2",
        f"{pdb2}.pdb",
        "-o",
        output_filename,
        "-m",
        "-ac",
        "-time",
        "-b",
    ]
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        output = result.stdout
        for line in output.splitlines():
            if "total time" in line:
                duration = float(line.split()[-1])
        return duration
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Command: {' '.join(command)}")
        print(f"Standard Output:\n{e.stdout}")
        print(f"Standard Error:\n{e.stderr}")


def extract_pdb(pair):
    pdb1 = pair["p_pdb"]
    pdb2 = pair["q_pdb"]
    pdb1_chain_id = pdb1.split("_")[-1]
    pdb2_chain_id = pdb2.split("_")[-1]
    pdb1_filename = f"pdb{pdb1.split('_')[0]}.ent.gz"
    pdb2_filename = f"pdb{pdb2.split('_')[0]}.ent.gz"
    print(pdb1_filename)
    if os.path.exists(pdb1_filename) and os.path.exists(pdb2_filename):
        pdb1_pdb = PandasPdb().read_pdb(pdb1_filename)
        pdb2_pdb = PandasPdb().read_pdb(pdb2_filename)
        pdb1_df = pdb1_pdb.df["ATOM"]
        pdb2_df = pdb2_pdb.df["ATOM"]
        pdb1_df = pdb1_df[pdb1_df["chain_id"] == pdb1_chain_id]
        pdb2_df = pdb2_df[pdb2_df["chain_id"] == pdb2_chain_id]
        pdb1_output_filename = f"{pdb1}.pdb"
        pdb2_output_filename = f"{pdb2}.pdb"
        pdb1_pdb.df["ATOM"] = pdb1_df
        pdb2_pdb.df["ATOM"] = pdb2_df
        print(pdb1_output_filename)
        pdb1_pdb.to_pdb(pdb1_output_filename)
        pdb2_pdb.to_pdb(pdb2_output_filename)


if __name__ == "__main__":
    main()
