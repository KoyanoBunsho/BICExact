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
        "fatcat_simulation_execution_time_improved.csv", index=False
    )


def run_fatcat(pdb1, pdb2, output_dir):
    output_filename = f"{output_dir}/{os.path.basename(pdb1)}_{os.path.basename(pdb2)}"
    command = [
        "./FATCAT_speed",
        "-p1", pdb1,
        "-p2", pdb2,
        "-o", output_filename,
        "-m",
        "-ac",
        "-time",
        "-b"
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
