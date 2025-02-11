import pandas as pd
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
        "dyndom_simulation_result/dyndom_simulation_execution_time.csv",
    )
    simulation_fatcat_df = read_fatcat_dyndom_data(
        "fatcat_simulation_hinge_count_result.csv",
        "fatcat_simulation_execution_time.csv",
    )
    df_dict = {}
    df_dict["BIC Exact"] = simulation_rmsdh_bic_df
    df_dict["BIC LH"] = simulation_fast_rmsdh_bic_df
    df_dict["AIC Exact"] = simulation_rmsdh_aic_df
    df_dict["AIC LH"] = simulation_fast_rmsdh_aic_df
    df_dict["DynDom"] = simulation_dyndom_df
    df_dict["FATCAT"] = simulation_fatcat_df
    df_dict["SE"] = simulation_auto_shibuya_df
    preprocess_and_merge(
        df_dict,
        key="DynDom",
        reference_key="BIC Exact",
        target_column="q_pdb",
        merge_columns=["actual_hinge_cnt_path", "actual_hinge_indices"],
    )
    preprocess_and_merge(
        df_dict,
        key="FATCAT",
        reference_key="BIC Exact",
        target_column="q_pdb",
        merge_columns=["actual_hinge_cnt_path", "actual_hinge_indices"],
    )
    plot_hinge_detection_accuracy(df_dict=df_dict)


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


def preprocess_and_merge(df_dict, key, reference_key, target_column, merge_columns):
    df_dict[key]["actual_hinge_cnt_path"] = (
        df_dict[key][target_column]
        .apply(lambda x: x.split("/")[1])
        .str.replace(".pdb", "", regex=False)
    )
    df_dict[key] = pd.merge(
        df_dict[key],
        df_dict[reference_key][merge_columns],
        on="actual_hinge_cnt_path",
        how="left",
    )


def plot_hinge_detection_accuracy(df_dict):
    plt.rcParams.update(
        {
            "font.size": 18,
            "axes.titlesize": 18,
            "axes.labelsize": 18,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
            "legend.fontsize": 12,
            "lines.markersize": 10,
        }
    )
    marker_styles = {
        "BIC Exact": "s",
        "BIC LH": "+",
        "FATCAT": "*",
        "DynDom": "v",
        "AIC Exact": "<",
        "AIC LH": ">",
        "SE": "X",
    }
    colors = plt.cm.get_cmap("tab10", 7)
    color_map = {
        "BIC Exact": colors(0),
        "BIC LH": colors(1),
        "FATCAT": colors(2),
        "DynDom": colors(3),
        "AIC Exact": colors(4),
        "AIC LH": colors(5),
        "SE": colors(6),
    }

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(24, 8))
    f_measure_label = "F-measure"
    precision_label = "Precision"
    recall_label = "Recall"
    f_measure_order = [
        "BIC Exact",
        "BIC LH",
        "FATCAT",
        "AIC LH",
        "AIC Exact",
        "DynDom",
        "SE",
    ]
    precision_order = [
        "BIC Exact",
        "BIC LH",
        "FATCAT",
        "DynDom",
        "AIC LH",
        "AIC Exact",
        "SE",
    ]
    # TODO: バグを修正する
    recall_order = [
        "SE",
        "AIC Exact",
        "AIC LH",
        "BIC Exact",
        "BIC LH",
        "FATCAT",
        "DynDom",
    ]
    f_measure_lines = []
    precision_lines = []
    recall_lines = []
    results = []
    results_precision = []
    results_recall = []
    all_f_measure = {}
    all_precision = {}
    all_recall = {}

    for _, (acc_method, precision_method, recall_method) in enumerate(
        zip(f_measure_order, precision_order, recall_order)
    ):
        precision_res = []
        recall_res = []
        f_measure_res = []
        acc_df = df_dict[acc_method]
        precision_df = df_dict[precision_method]
        recall_df = df_dict[recall_method]
        res_acc = calc_acc_df(
            acc_df[acc_df["actual_hinge_cnt"] != 0]["actual_hinge_indices"]
            .reset_index()
            .rename(columns={"actual_hinge_indices": "hinge_index"}),
            acc_df[acc_df["actual_hinge_cnt"] != 0]["hinge_index"].reset_index(),
        )
        res_precision = calc_acc_df(
            precision_df[precision_df["actual_hinge_cnt"] != 0]["actual_hinge_indices"]
            .reset_index()
            .rename(columns={"actual_hinge_indices": "hinge_index"}),
            precision_df[precision_df["actual_hinge_cnt"] != 0][
                "hinge_index"
            ].reset_index(),
        )
        res_recall = calc_acc_df(
            recall_df[recall_df["actual_hinge_cnt"] != 0]["actual_hinge_indices"]
            .reset_index()
            .rename(columns={"actual_hinge_indices": "hinge_index"}),
            recall_df[recall_df["actual_hinge_cnt"] != 0]["hinge_index"].reset_index(),
        )
        f_measure = res_acc["F-measure_distance_3"]
        precision = res_precision["Precision_distance_3"]
        recall = res_recall["Recall_distance_3"]
        all_f_measure[acc_method] = f_measure
        all_precision[precision_method] = precision
        all_recall[recall_method] = recall

        for cnt in range(1, 6):
            tmp_acc = acc_df[acc_df["actual_hinge_cnt"] == cnt]
            tmp_precision = precision_df[precision_df["actual_hinge_cnt"] == cnt]
            tmp_recall = recall_df[recall_df["actual_hinge_cnt"] == cnt]
            res_acc = calc_acc_df(
                tmp_acc["actual_hinge_indices"]
                .reset_index()
                .rename(columns={"actual_hinge_indices": "hinge_index"}),
                tmp_acc["hinge_index"].reset_index(),
            )
            res_precision = calc_acc_df(
                tmp_precision["actual_hinge_indices"]
                .reset_index()
                .rename(columns={"actual_hinge_indices": "hinge_index"}),
                tmp_precision["hinge_index"].reset_index(),
            )
            res_recall = calc_acc_df(
                tmp_recall["actual_hinge_indices"]
                .reset_index()
                .rename(columns={"actual_hinge_indices": "hinge_index"}),
                tmp_recall["hinge_index"].reset_index(),
            )
            f_measure = res_acc["F-measure_distance_3"]
            precision = res_precision["Precision_distance_3"]
            recall = res_recall["Recall_distance_3"]
            f_measure_res.append(f_measure)
            precision_res.append(precision)
            recall_res.append(recall)
            results.append(
                {
                    "Method": acc_method,
                    "k": cnt,
                    "F-measure": f_measure,
                }
            )
            results_precision.append(
                {
                    "Method": precision_method,
                    "k": cnt,
                    "Precision": precision,
                }
            )
            results_recall.append(
                {
                    "Method": recall_method,
                    "k": cnt,
                    "Recall": recall,
                }
            )
        (f_measure_line,) = ax1.plot(
            range(1, 6),
            f_measure_res,
            label=acc_method,
            marker=marker_styles[acc_method],
            color=color_map[acc_method],
        )
        (precision_line,) = ax2.plot(
            range(1, 6),
            precision_res,
            label=precision_method,
            marker=marker_styles[precision_method],
            color=color_map[precision_method],
        )
        (recall_line,) = ax3.plot(
            range(1, 6),
            recall_res,
            label=recall_method,
            marker=marker_styles[recall_method],
            color=color_map[recall_method],
        )
        f_measure_lines.append(f_measure_line)
        precision_lines.append(precision_line)
        recall_lines.append(recall_line)

    ax1.text(
        0.5,
        -0.2,
        f"({chr(97 + 0)}) F-measure of all methods against the true hinge number",
        transform=ax1.transAxes,
        fontsize=18,
        va="center",
        ha="center",
    )
    ax2.text(
        0.5,
        -0.2,
        f"({chr(97 + 1)}) Precision of all methods against the true hinge number",
        transform=ax2.transAxes,
        fontsize=18,
        va="center",
        ha="center",
    )
    ax3.text(
        0.5,
        -0.2,
        f"({chr(97 + 2)}) Recall of all methods against the true hinge number",
        transform=ax3.transAxes,
        fontsize=18,
        va="center",
        ha="center",
    )
    x_ticks = range(1, 6)
    ax1.set_xlim(1, 5)
    ax1.set_xticks(x_ticks)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("True hinge number")
    ax1.set_ylabel(f_measure_label)
    ax1.legend(f_measure_lines, f_measure_order, loc="best")
    ax2.set_xlim(1, 5)
    ax2.set_xticks(x_ticks)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel(precision_label)
    ax2.set_xlabel("True hinge number")
    ax2.legend(precision_lines, precision_order, loc="best")
    ax3.set_xlim(1, 5)
    ax3.set_xticks(x_ticks)
    ax3.set_ylim(0, 1)
    ax3.set_ylabel(recall_label)
    ax3.set_xlabel("True hinge number")
    ax3.legend(recall_lines, recall_order, loc="best")
    fig.tight_layout()
    plt.savefig(f"figures/{f_measure_label}_simulation_data.png")
    plt.savefig(f"figures/{f_measure_label}_simulation_data.svg", format="svg")
    plt.clf()
    plt.close()
    results_df = pd.DataFrame(results)
    results_df.to_csv("figures/f_measure_results.csv", index=False)
    results_precision_df = pd.DataFrame(results_precision)
    results_precision_df.to_csv("figures/precision_results.csv", index=False)
    results_recall_df = pd.DataFrame(results_recall)
    results_recall_df.to_csv("figures/recall_results.csv", index=False)
    all_f_measure = dict(
        sorted(all_f_measure.items(), key=lambda item: item[1], reverse=True)
    )
    all_precision = dict(
        sorted(all_precision.items(), key=lambda item: item[1], reverse=True)
    )
    all_recall = dict(
        sorted(all_recall.items(), key=lambda item: item[1], reverse=True)
    )
    pd.DataFrame(all_f_measure, index=[0]).to_csv(
        "figures/all_f_measure_results.csv", index=False
    )
    pd.DataFrame(all_precision, index=[0]).to_csv(
        "figures/all_precision_results.csv", index=False
    )
    pd.DataFrame(all_recall, index=[0]).to_csv(
        "figures/all_recall_results.csv", index=False
    )


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


if __name__ == "__main__":
    main()
