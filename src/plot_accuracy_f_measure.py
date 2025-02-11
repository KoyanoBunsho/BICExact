import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"


def main():
    os.makedirs("figures", exist_ok=True)
    simulation_rmsdh_bic_df = read_simulation_data(
        "rmsdh_result/simulation_parallel_rmsdh_bic_0.894893output_simulation_file.csv"
    )
    simulation_fast_rmsdh_bic_df = read_simulation_data(
        "rmsdh_result/simulation_parallel_fast_rmsdh_bic_0.894893output_simulation_file.csv"
    )
    simulation_rmsdh_aic_df = read_simulation_data(
        "rmsdh_result/simulation_parallel_rmsdh_aic_0.894893output_simulation_file.csv"
    )
    simulation_fast_rmsdh_aic_df = read_simulation_data(
        "rmsdh_result/simulation_parallel_fast_rmsdh_aic_0.894893output_simulation_file.csv"
    )
    simulation_auto_shibuya_df = read_simulation_data(
        "rmsdh_result/simulation_auto_shibuya_parallel.csv"
    )
    simulation_dyndom_df = read_fatcat_dyndom_data(
        "dyndom_simulation_hinge_count_result.csv",
        "dyndom_simulation_result/dyndom_simulation_execution_time_improved.csv",
    )
    simulation_fatcat_df = read_fatcat_dyndom_data(
        "fatcat_simulation_hinge_count_result.csv",
        "fatcat_simulation_execution_time_improved.csv",
    )
    df_dict = {}
    df_dict["DP based on BIC"] = simulation_rmsdh_bic_df
    df_dict["LH DP based on BIC"] = simulation_fast_rmsdh_bic_df
    df_dict["DP based on AIC"] = simulation_rmsdh_aic_df
    df_dict["LH DP based on AIC"] = simulation_fast_rmsdh_aic_df
    df_dict["DynDom"] = simulation_dyndom_df
    df_dict["FATCAT"] = simulation_fatcat_df
    df_dict["SE"] = simulation_auto_shibuya_df
    df_dict["DynDom"]["actual_hinge_cnt_path"] = (
        df_dict["DynDom"]["q_pdb"]
        .apply(lambda x: x.split("/")[1])
        .str.replace(".pdb", "")
    )
    df_dict["FATCAT"]["actual_hinge_cnt_path"] = (
        df_dict["FATCAT"]["q_pdb"]
        .apply(lambda x: x.split("/")[1])
        .str.replace(".pdb", "")
    )
    df_dict["DynDom"] = pd.merge(
        df_dict["DynDom"],
        df_dict["DP based on BIC"][["actual_hinge_cnt_path", "actual_hinge_indices"]],
        on="actual_hinge_cnt_path",
        how="left",
    )
    df_dict["FATCAT"] = pd.merge(
        df_dict["FATCAT"],
        df_dict["DP based on BIC"][["actual_hinge_cnt_path", "actual_hinge_indices"]],
        on="actual_hinge_cnt_path",
        how="left",
    )
    for how in ["hinge_count", "hinge_detection"]:
        plot_hinge_detection_accuracy(df_dict=df_dict, how=how)
    measure_computation_time(df_dict=df_dict)


def read_simulation_data(file_path):
    df = (
        pd.read_csv(file_path)
        .fillna("")
        .rename(columns={"actual_hinge_cnt": "actual_hinge_cnt_path"})
    )
    df["actual_hinge_cnt_path"] = (
        df["actual_hinge_cnt_path"]
        .apply(lambda x: str(x).split("/")[1])
        .str.replace(".csv", "")
    )
    df["actual_hinge_cnt"] = (
        df["actual_hinge_cnt_path"]
        .str.split("_")
        .apply(lambda x: str(x[3]))
        .astype(int)
    )
    return df


def read_fatcat_dyndom_data(file_path1, file_path2):
    df1 = pd.read_csv(file_path1).fillna("")
    df1["actual_hinge_cnt_path"] = (
        df1["q_pdb"].apply(lambda x: str(x).split("/")[1]).str.replace(".pdb", "")
    )
    df1["actual_hinge_cnt"] = (
        df1["q_pdb"].str.split("_").apply(lambda x: str(x[4])).astype(int)
    )
    df1 = df1.rename(columns={"detected_hinge_count": "hinge_cnt"})
    df2 = pd.read_csv(file_path2)
    df1["exec_time (s)"] = df2["execution_time"]
    return df1


def plot_hinge_detection_accuracy(df_dict, how):
    plt.rcParams.update(
        {
            "font.size": 30,
            "axes.titlesize": 30,
            "axes.labelsize": 30,
            "xtick.labelsize": 30,
            "ytick.labelsize": 30,
            "legend.fontsize": 30,
            "lines.markersize": 10,
        }
    )
    marker_styles = {
        "DP based on BIC": "s",
        "LH DP based on BIC": "^",
        "FATCAT": "*",
        "DynDom": "v",
        "DP based on AIC": "<",
        "LH DP based on AIC": ">",
        "SE": "X",
    }
    colors = plt.cm.get_cmap("tab10", 7)
    color_map = {
        "DP based on BIC": colors(0),
        "LH DP based on BIC": colors(1),
        "FATCAT": colors(2),
        "DynDom": colors(3),
        "DP based on AIC": colors(4),
        "LH DP based on AIC": colors(5),
        "SE": colors(6),
    }

    # 3列レイアウトに変更（凡例用の幅を調整）
    fig, axes = plt.subplots(
        nrows=1, ncols=3, figsize=(18, 7), gridspec_kw={"width_ratios": [1, 1, 0.2]}
    )
    ax1, ax2, ax3 = axes

    if how == "hinge_position":
        accuracy_label = "F-measure"
        mae_label = "Precision"
    else:
        accuracy_label = "Accuracy"
        mae_label = "MAE"

    accuracy_order = [
        "DP based on BIC",
        "LH DP based on BIC",
        "DP based on AIC",
        "LH DP based on AIC",
        "FATCAT",
        "DynDom",
        "SE",
    ]

    for method in accuracy_order:
        acc_df = df_dict[method]
        mae_df = df_dict[method]

        accuracy_res = []
        mae_res = []

        for cnt in range(0, 6):
            tmp_acc = acc_df[acc_df["actual_hinge_cnt"] == cnt]
            tmp_mae = mae_df[mae_df["actual_hinge_cnt"] == cnt]

            accuracy = (
                len(tmp_acc[tmp_acc["actual_hinge_cnt"] == tmp_acc["hinge_cnt"]])
                / len(tmp_acc)
                if len(tmp_acc) > 0
                else 0
            )
            mae = (
                mean_absolute_error(tmp_mae["actual_hinge_cnt"], tmp_mae["hinge_cnt"])
                if len(tmp_mae) > 0
                else 0
            )

            accuracy_res.append(accuracy)
            mae_res.append(mae)

        # Accuracy plot
        ax1.plot(
            range(0, 6),
            accuracy_res,
            label=method,
            marker=marker_styles[method],
            color=color_map[method],
        )

        # MAE plot
        ax2.plot(
            range(0, 6),
            mae_res,
            label=method,
            marker=marker_styles[method],
            color=color_map[method],
        )

    ax1.set_xlabel("True hinge number")
    ax1.set_ylabel(accuracy_label)
    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 1)
    ax2.set_xlim(0, 5)
    ax2.set_ylim(0)
    ax2.set_xlabel("True hinge number")
    ax2.set_ylabel(mae_label)
    ax1.text(
        0.5,
        -0.3,
        f"({chr(97 + 0)})",
        transform=ax1.transAxes,
        fontsize=30,
        va="center",
        ha="center",
    )
    ax2.text(
        0.5,
        -0.3,
        f"({chr(97 + 1)})",
        transform=ax2.transAxes,
        fontsize=30,
        va="center",
        ha="center",
    )
    # 3列目（凡例専用）の設定
    ax3.axis("off")
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker=marker_styles[method],
            color="w",
            markerfacecolor=color_map[method],
            markersize=10,
            label=method,
        )
        for method in accuracy_order
    ]
    ax3.legend(handles=handles, loc="center", frameon=False, title="Methods")
    plt.tight_layout()
    plt.savefig(
        f"figures/{accuracy_label}_simulation_data_{how}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        f"figures/{accuracy_label}_simulation_data_{how}.svg",
        format="svg",
        bbox_inches="tight",
    )
    plt.close()


def plot_hinge_detection_accuracy_old(df_dict, how):
    plt.rcParams.update(
        {
            "font.size": 30,
            "axes.titlesize": 30,
            "axes.labelsize": 30,
            "xtick.labelsize": 30,
            "ytick.labelsize": 30,
            "legend.fontsize": 30,
            "lines.markersize": 10,
        }
    )
    marker_styles = {
        "DP based on BIC": "s",
        "LH DP based on BIC": "^",
        "FATCAT": "*",
        "DynDom": "v",
        "DP based on AIC": "<",
        "LH DP based on AIC": ">",
        "SE": "X",
    }
    colors = plt.cm.get_cmap("tab10", 7)
    color_map = {
        "DP based on BIC": colors(0),
        "LH DP based on BIC": colors(1),
        "FATCAT": colors(2),
        "DynDom": colors(3),
        "DP based on AIC": colors(4),
        "LH DP based on AIC": colors(5),
        "SE": colors(6),
    }

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(30, 6))
    if how == "hinge_position":
        accuracy_label = "F-measure"
        mae_label = "Precision"
    else:
        accuracy_label = "Accuracy"
        mae_label = "MAE"

    accuracy_order = [
        "DP based on BIC",
        "LH DP based on BIC",
        "FATCAT",
        "DynDom",
        "DP based on AIC",
        "LH DP based on AIC",
        "SE",
    ]
    mae_order = [
        "SE",
        "DP based on AIC",
        "LH DP based on AIC",
        "DynDom",
        "FATCAT",
        "DP based on BIC",
        "LH DP based on BIC",
    ]

    accuracy_lines = []
    mae_lines = []
    results = []
    results_mae = []

    for _, (acc_method, mae_method) in enumerate(zip(accuracy_order, mae_order)):
        mae_res = []
        rmse_res = []
        accuracy_res = []
        acc_df = df_dict[acc_method]
        mae_df = df_dict[mae_method]
        if how == "hinge_count":
            accuracy = (
                len(acc_df[acc_df["actual_hinge_cnt"] == acc_df["hinge_cnt"]])
                / len(acc_df)
                if len(acc_df) > 0
                else 0
            )
            mae = (
                mean_absolute_error(mae_df["actual_hinge_cnt"], mae_df["hinge_cnt"])
                if len(mae_df) > 0
                else 0
            )
            print(f"{acc_method}: Overall Accuracy, {format_stats(accuracy)}")
            print(f"{mae_method}: Overall MAE, {format_stats(mae)}")
        for cnt in range(0, 6):
            tmp_acc = acc_df[acc_df["actual_hinge_cnt"] == cnt]
            tmp_mae = mae_df[mae_df["actual_hinge_cnt"] == cnt]
            if how == "hinge_count":
                accuracy = (
                    len(tmp_acc[tmp_acc["actual_hinge_cnt"] == tmp_acc["hinge_cnt"]])
                    / len(tmp_acc)
                    if len(tmp_acc) > 0
                    else 0
                )
                mae = (
                    mean_absolute_error(
                        tmp_mae["actual_hinge_cnt"], tmp_mae["hinge_cnt"]
                    )
                    if len(tmp_mae) > 0
                    else 0
                )
                rmse = (
                    mean_squared_error(
                        tmp_acc["actual_hinge_cnt"], tmp_acc["hinge_cnt"]
                    )
                    ** 0.5
                    if len(tmp_acc) > 0
                    else 0
                )
            else:
                res_acc = calc_acc_df(
                    tmp_acc["actual_hinge_indices"]
                    .reset_index()
                    .rename(columns={"actual_hinge_indices": "hinge_index"}),
                    tmp_acc["hinge_index"].reset_index(),
                )
                res_mae = calc_acc_df(
                    tmp_mae["actual_hinge_indices"]
                    .reset_index()
                    .rename(columns={"actual_hinge_indices": "hinge_index"}),
                    tmp_mae["hinge_index"].reset_index(),
                )
                accuracy = res_acc["F-measure_distance_3"]
                mae = res_mae["Precision_distance_3"]
                rmse = res_acc["Recall_distance_3"]
            accuracy_res.append(accuracy)
            mae_res.append(mae)
            rmse_res.append(rmse)
            results.append(
                {
                    "Method": acc_method,
                    "k": cnt,
                    "Accuracy": accuracy,
                }
            )
            results_mae.append(
                {
                    "Method": mae_method,
                    "k": cnt,
                    "Accuracy": mae,
                }
            )

        (accuracy_line,) = ax1.plot(
            range(0, 6),
            accuracy_res,
            label=acc_method,
            marker=marker_styles[acc_method],
            color=color_map[acc_method],
        )
        (mae_line,) = ax2.plot(
            range(0, 6),
            mae_res,
            label=mae_method,
            marker=marker_styles[mae_method],
            color=color_map[mae_method],
        )
        accuracy_lines.append(accuracy_line)
        mae_lines.append(mae_line)

    ax1.text(
        0.5,
        -1,
        f"({chr(97 + 0)}) Accuracy against the true hinge number",
        transform=ax1.transAxes,
        fontsize=30,
        va="center",
        ha="center",
    )
    ax2.text(
        0.5,
        -1,
        f"({chr(97 + 1)}) MAE against the true hinge number",
        transform=ax2.transAxes,
        fontsize=30,
        va="center",
        ha="center",
    )
    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("True hinge number")
    ax1.set_ylabel(accuracy_label)
    ax1.legend(accuracy_lines, accuracy_order, loc="best")
    ax2.set_xlim(0, 5)
    if how == "hinge_position":
        ax2.set_ylim(0, 1)
    else:
        ax2.set_ylim(0)
    ax2.set_ylabel(mae_label)
    ax2.set_xlabel("True hinge number")
    ax2.legend(mae_lines, mae_order, loc="best")
    fig.tight_layout()
    plt.savefig(f"figures/{accuracy_label}_simulation_data_{how}.png")
    plt.savefig(f"figures/{accuracy_label}_simulation_data_{how}.svg", format="svg")
    plt.clf()
    plt.close()
    results_df = pd.DataFrame(results)
    results_df.to_csv(f"figures/accuracy_results_{how}.csv", index=False)
    results_mae_df = pd.DataFrame(results_mae)
    results_mae_df.to_csv(f"figures/mae_results_{how}.csv", index=False)


def calc_ans_dyndom(exp, detect, d=0):
    true_hinge_indices = exp.split(" : ")
    detected_hinge_indices = detect.split(" : ")
    TP = 0
    FP = 0
    FN = 0
    if detected_hinge_indices != [""]:
        detected_ranges = [
            (int(label) - d, int(label) + d) for label in detected_hinge_indices
        ]
    else:
        detected_ranges = []
    for true in true_hinge_indices:
        if true != "":
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        else:
            if detected_hinge_indices == [""]:
                TP += 1
    for lower, upper in detected_ranges:
        if true != "":
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        else:
            FP += 1
    return {"TP": TP, "FP": FP, "FN": FN}


def calc_acc_df(df, heuristic_df):
    df["hinge_index"] = df["hinge_index"].fillna("")
    heuristic_df["hinge_index"] = heuristic_df["hinge_index"].fillna("")
    f_measure_dict = {}
    # バグがないか検算する
    # f-measure, precision, recallを全て載せる
    for d in [3]:
        acc = []
        for i in range(len(df)):
            acc.append(
                calc_ans_dyndom(
                    df.loc[i]["hinge_index"], heuristic_df.loc[i]["hinge_index"], d
                )
            )
        acc_df = pd.DataFrame(acc)
        TP = acc_df["TP"].sum()
        FP = acc_df["FP"].sum()
        FN = acc_df["FN"].sum()
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        f_measure_dict[f"F-measure_distance_{d}"] = f_measure
        f_measure_dict[f"Precision_distance_{d}"] = precision
        f_measure_dict[f"Recall_distance_{d}"] = recall
        f_measure_dict[f"TP_{d}"] = TP
        f_measure_dict[f"FP_{d}"] = FP
        f_measure_dict[f"FN_{d}"] = FN
    return f_measure_dict


def format_stats(stats):
    return f"{stats:.3f}"


def measure_computation_time(df_dict):
    for method, df in df_dict.items():
        print(
            f"{method}: total computation time is {format_stats(df['exec_time (s)'].sum())}"
        )
        print(
            f"{method}: average computation time is {format_stats(df['exec_time (s)'].mean())}"
        )


if __name__ == "__main__":
    main()
