import os
import pandas as pd
import subprocess
import time
from joblib import Parallel, delayed
from tqdm import tqdm
import glob


def main():
    result_directory = "dyndom_simulation_result"
    os.makedirs(result_directory, exist_ok=True)
    simulation_file = "output_simulation_file_for_fatcat.csv"
    df = pd.read_csv(simulation_file)
    make_dyndom_command_files(df, result_directory)
    durations = run_all_dyndom(result_directory)
    execution_time_list = [
        {"pdb_pair": os.path.basename(pdb_pair), "execution_time": duration}
        for pdb_pair, duration in durations
    ]
    pd.DataFrame(execution_time_list).to_csv(
        os.path.join(result_directory, "dyndom_simulation_execution_time_improved.csv"),
        index=False,
    )


def make_dyndom_command_files(df, result_directory):
    print("Start making dyndom command files")
    Parallel(n_jobs=os.cpu_count())(
        delayed(make_dyndom_command_file)(row["p_pdb"], row["q_pdb"], result_directory)
        for _, row in tqdm(df.iterrows(), total=len(df))
    )


def make_dyndom_command_file(pdb1, pdb2, result_directory):
    pdb1_basename = os.path.basename(pdb1).replace(".pdb", "")
    pdb2_basename = os.path.basename(pdb2).replace(".pdb", "")
    title = f"{pdb1_basename}_{pdb2_basename}"
    command_file_path = os.path.join(result_directory, f"{title}.command")
    chain1id = pdb1_basename.split("_")[1]
    chain2id = pdb2_basename.split("_")[1]
    command_file_content = (
        f"title={title}\n"
        f"filename1=./simulation_data/{os.path.basename(pdb1)}\n"
        f"chain1id={chain1id}\n"
        f"filename2=./simulation_data/{os.path.basename(pdb2)}\n"
        f"chain2id={chain2id}\n"
        "window=5\n"
        "domain=20\n"
        "ratio=1.0\n"
    )
    with open(command_file_path, "w") as file:
        file.write(command_file_content)


def run_dyndom(command_file_path):
    command = ["./DynDom_speed", f"{command_file_path}"]
    start_time = time.time()
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        output = result.stdout
        duration = None
        for line in output.splitlines():
            if "Time for CLUSTER2:" in line:
                duration = float(line.split(": ")[1])
                print(duration)
        if duration is None:
            duration = time.time() - start_time
        return duration
    except subprocess.CalledProcessError as e:
        print(f"Error in executing command: {command}\nError: {str(e)}")
        return time.time() - start_time


def run_all_dyndom(result_directory):
    command_files = glob.glob(f"{result_directory}/*.command")
    results = Parallel(n_jobs=os.cpu_count())(
        delayed(run_dyndom_with_print)(command_file)
        for command_file in tqdm(command_files, total=len(command_files))
    )
    return results


def run_dyndom_with_print(command_file):
    duration = run_dyndom(command_file)
    return (command_file, duration)


if __name__ == "__main__":
    main()
