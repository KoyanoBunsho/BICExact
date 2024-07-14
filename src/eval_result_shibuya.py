import pandas as pd
import os
from sklearn.metrics import mean_absolute_error, accuracy_score

expert_annotated_list = [
    "35",
    "25 : 55",
    "98 : 110",
    "13 : 81",
    "89 : 264",
    "135",
    "150 : 318",
    "86 : 177",
    "91 : 251",
    "91 : 193",
    "104 : 236",
    "36 : 63 : 102",
]

fatcat_list = ["", "", "", "", "91", "", "", "89 : 183", "92 : 252", "81", "99", ""]
dyndom_list = [
    "",
    "32 : 35 : 42 : 59 : 74 : 84",  # 31 - 32, 34 - 35, 41 - 42, 58 - 59, 73 - 74, 82 - 85
    "96 : 119",  # 94 - 97, 114 - 123
    "11 : 75",  # 9 - 12, 71 - 78
    "180 : 186 : 265",  # 179 - 180, 185 - 186, 264 - 265
    "83 : 100 : 117 : 135 : 187 : 196",  # 82 - 83, 99 - 100, 113 - 121, 134 - 135, 186 - 188, 195 - 197
    "153 : 171 : 205 : 218 : 230 : 243 : 318 : 343 : 368 : 371",  # 152 - 153, 170 - 171, 203 - 206, 214 - 221, 229 - 230, 242 - 243, 317 - 319, 342 - 343, 367 - 368, 370 - 371
    "88 : 183",  # 87 - 89, 181 - 184
    "91 : 251 : 330 : 683",  # 89 - 92, 249 - 253, 329 - 330, 682 - 683
    "91 : 192",  # 89 - 92, 190 - 193
    "103 : 235 : 265",  # 101 - 104, 233 - 237, 263 - 267
    "40 : 70",  # 38 - 41, 66 - 73
]


def main():
    os.makedirs("figures", exist_ok=True)
    auto_shibuya = pd.read_csv("rmsdh_result/shibuya_auto_rmsdhk.csv").fillna("")
    rmsdh_bic = pd.read_csv("rmsdh_result/shibuya_rmsdh_bic_0.894893.csv").fillna("")
    rmsdh_aic = pd.read_csv("rmsdh_result/shibuya_rmsdh_aic_0.894893.csv").fillna("")
    fast_rmsdh_bic = pd.read_csv(
        "rmsdh_result/shibuya_fast_rmsdh_bic_0.894893.csv"
    ).fillna("")
    fast_rmsdh_aic = pd.read_csv(
        "rmsdh_result/shibuya_fast_rmsdh_aic_0.894893.csv"
    ).fillna("")
    dyndom_execution_time = pd.read_csv(
        "all_pdb/dyndom_execution_time_shibuya.csv"
    ).rename(columns={"execution_time": "exec_time (s)"})
    fatcat_execution_time = pd.read_csv(
        "all_pdb/fatcat_execution_time_shibuya.csv"
    ).rename(columns={"execution_time": "exec_time (s)"})
    methods = {
        "aic": rmsdh_aic["hinge_index"].tolist(),
        "lh_aic": fast_rmsdh_aic["hinge_index"].tolist(),
        "bic": rmsdh_bic["hinge_index"].tolist(),
        "lh_bic": fast_rmsdh_bic["hinge_index"].tolist(),
        "SE": auto_shibuya["hinge_index"].to_list(),
        "FATCAT": fatcat_list,
        "DynDom": dyndom_list,
    }
    comp_time_dict = {
        "BIC Exact": rmsdh_bic,
        "BIC LH": fast_rmsdh_bic,
        "SE": auto_shibuya,
        "FATCAT": fatcat_execution_time,
        "DynDom": dyndom_execution_time,
    }
    results = {}
    for name, method_list in methods.items():
        results[name] = calc_stats(method_list, expert_annotated_list)
    print(
        f"BIC: {format_stats(results['bic'])}\n"
        f"BIC LH: {format_stats(results['lh_bic'])}\n"
        f"AIC: {format_stats(results['aic'])}\n"
        f"AIC LH: {format_stats(results['lh_aic'])}\n"
        f"SE: {format_stats(results['SE'])}\n"
        f"FATCAT: {format_stats(results['FATCAT'])}\n"
        f"DynDom: {format_stats(results['DynDom'])}\n"
    )
    measure_comp_time(comp_time_dict)


def format_stats(stats, precision=2):
    return f"({stats[0]:.{precision}f}, {stats[1]:.{precision}f})"


def format_precision(stats, precision=3):
    return f"{stats:.{precision}f}"


def calc_stats(method_list, expert_list):
    method_hinge_counts = [len(hinges.split(" : ")) for hinges in method_list]
    expert_hinge_counts = [len(hinge.split(" : ")) for hinge in expert_list]
    mae = mean_absolute_error(expert_hinge_counts, method_hinge_counts)
    accuracy = accuracy_score(expert_hinge_counts, method_hinge_counts)
    return (accuracy, mae)


def measure_comp_time(comp_time_dict):
    for method, df in comp_time_dict.items():
        print(
            f"{method}: total computation time is {format_precision(df['exec_time (s)'].sum())}"
        )
        print(
            f"{method}: average computation time is {format_precision(df['exec_time (s)'].mean())}"
        )


if __name__ == "__main__":
    main()
