import os
import subprocess
import time
import pandas as pd
from joblib import Parallel, delayed


def main():
    output_dir = "fatcat_simulation_result"
    os.makedirs(output_dir, exist_ok=True)
    simulation_file = "output_simulation_file_for_fatcat.csv"
    df = pd.read_csv(simulation_file)
    durations = compare_pdb_files(df, output_dir)
    execution_time_list = []
    for pdb, duration in durations:
        execution_time_list.append({"pdb": pdb, "execution_time": duration})
    pd.DataFrame(execution_time_list).to_csv(
        "fatcat_simulation_execution_time.csv", index=False
    )


def run_fatcat(pdb1, pdb2, output_dir):
    output_filename = f"{output_dir}/{os.path.basename(pdb1)}_{os.path.basename(pdb2)}"
    command = f"./FATCAT -p1 {pdb1} -p2 {pdb2} -o {output_filename} -m -ac -t"
    start_time = time.time()
    try:
        subprocess.run(command, shell=True, check=True)
        duration = time.time() - start_time
        print(f"Command executed successfully: {command}")
        return duration
    except subprocess.CalledProcessError as e:
        print(f"Error in executing command: {command}\nError: {str(e)}")
        return time.time() - start_time


def compare_pdb_files(df, output_dir):
    durations = Parallel(n_jobs=os.cpu_count())(
        delayed(process_pair)(row["p_pdb"], row["q_pdb"], output_dir)
        for _, row in df.iterrows()
    )
    return durations


def process_pair(pdb1, pdb2, output_dir):
    duration = run_fatcat(pdb1, pdb2, output_dir)
    return (os.path.basename(pdb2), duration)


if __name__ == "__main__":
    main()
