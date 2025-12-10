#!/usr/bin/env python3
"""
Analyze benchmarking results for Linear Probing vs Robin Hood hashing.

- Reads results.csv (with capacity_log2, max_probes) from the current directory.
- Optionally reads probe_hist.csv for probe-length distributions.
- Computes time_per_op.
- Aggregates by (capacity_log2, table_type, load_factor, phase).
- Generates plots:
  * Time per op vs load factor (all capacities on one plot per phase).
  * Avg probes vs load factor (all capacities on one plot per phase).
  * Max probes vs load factor (all capacities on one plot per phase).
  * Time vs probes (all capacities on one plot per phase).
  * Boxplot of time_per_op across repetitions (lookup_fail), faceted by capacity.
  * Probe quantiles (p50/p90/p99) vs load factor (from probe_hist.csv) per phase.
- Saves PNGs in the current directory.
- Prints aggregated summary to console.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 140,
    "axes.grid": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})
sns.set_style("whitegrid")

PALETTE = {"linear": "#1f77b4", "robinhood": "#d62728"}
TABLE_ORDER = ["linear", "robinhood"]


def load_results(csv_path: Path) -> pd.DataFrame:
    # Try different encodings
    encodings = ['utf-8', 'utf-8-sig', 'utf-16', 'latin-1']
    for enc in encodings:
        try:
            df = pd.read_csv(csv_path, encoding=enc)
            df["table_type"] = df["table_type"].astype(str).str.strip().str.lower()
            df["phase"] = df["phase"].astype(str).str.strip()
            df["load_factor"] = df["load_factor"].astype(float)
            df["time_per_op"] = df["time_seconds"] / df["num_ops"]
            return df
        except (UnicodeDecodeError, UnicodeError):
            continue
    raise ValueError(f"Could not read {csv_path} with any supported encoding")


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    agg = (
        df.groupby(["capacity_log2", "table_type", "load_factor", "phase"])
        .agg(
            time_per_op_mean=("time_per_op", "mean"),
            time_per_op_std=("time_per_op", "std"),
            avg_probes_mean=("avg_probes", "mean"),
            avg_probes_std=("avg_probes", "std"),
            max_probes_mean=("max_probes", "mean"),
            max_probes_std=("max_probes", "std"),
            reps=("time_per_op", "count"),
        )
        .reset_index()
        .sort_values(["phase", "capacity_log2", "load_factor", "table_type"])
    )
    agg["time_per_op_std"] = agg["time_per_op_std"].fillna(0.0)
    agg["avg_probes_std"] = agg["avg_probes_std"].fillna(0.0)
    agg["max_probes_std"] = agg["max_probes_std"].fillna(0.0)
    return agg


def plot_metric_vs_load(agg: pd.DataFrame, metric_mean: str, metric_std: str,
                        ylabel: str, prefix: str, title_prefix: str) -> None:
    phases = agg["phase"].unique()
    capacities = sorted(agg["capacity_log2"].unique())

    for phase in phases:
        phase_df = agg[agg["phase"] == phase]
        fig, ax = plt.subplots(figsize=(7, 4))
        for cap in capacities:
            cap_df = phase_df[phase_df["capacity_log2"] == cap]
            for tbl in TABLE_ORDER:
                sub = cap_df[cap_df["table_type"] == tbl].sort_values("load_factor")
                if sub.empty:
                    continue
                label = f"{tbl} (2^{cap})"
                ax.errorbar(
                    sub["load_factor"],
                    sub[metric_mean],
                    yerr=sub[metric_std],
                    label=label,
                    marker="o",
                    capsize=3,
                    color=PALETTE.get(tbl, None),
                    linestyle="--" if cap % 2 else "-",
                )
        ax.set_xlabel("Load factor")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_prefix}: {phase}")
        ax.legend(title="Table / capacity", ncol=2)
        fig.tight_layout()
        fig.savefig(f"{prefix}_{phase}.png")
        plt.close(fig)


def plot_time_vs_probes(agg: pd.DataFrame) -> None:
    phases = agg["phase"].unique()
    capacities = sorted(agg["capacity_log2"].unique())
    for phase in phases:
        phase_df = agg[agg["phase"] == phase]
        fig, ax = plt.subplots(figsize=(7, 4))
        for cap in capacities:
            cap_df = phase_df[phase_df["capacity_log2"] == cap]
            for tbl in TABLE_ORDER:
                sub = cap_df[cap_df["table_type"] == tbl].sort_values("avg_probes_mean")
                if sub.empty:
                    continue
                label = f"{tbl} (2^{cap})"
                ax.plot(
                    sub["avg_probes_mean"],
                    sub["time_per_op_mean"],
                    label=label,
                    marker="o",
                    color=PALETTE.get(tbl, None),
                    linestyle="--" if cap % 2 else "-",
                )
        ax.set_xlabel("Average probes per operation")
        ax.set_ylabel("Time per operation (s)")
        ax.set_title(f"Time vs Probes: {phase}")
        ax.legend(title="Table / capacity", ncol=2)
        fig.tight_layout()
        fig.savefig(f"time_vs_probes_{phase}.png")
        plt.close(fig)


def plot_boxplot_lookup_fail(df: pd.DataFrame) -> None:
    phase = "lookup_fail"
    sub = df[df["phase"] == phase]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.boxplot(
        data=sub,
        x="load_factor",
        y="time_per_op",
        hue="table_type",
        palette=PALETTE,
        ax=ax,
    )
    ax.set_xlabel("Load factor")
    ax.set_ylabel("Time per operation (s)")
    ax.set_title("Time per operation distribution (lookup_fail)")
    ax.legend(title="Table type")
    fig.tight_layout()
    fig.savefig("boxplot_time_per_op_lookup_fail.png")
    plt.close(fig)


def load_probe_hist(csv_path: Path) -> pd.DataFrame | None:
    if not csv_path.exists():
        return None
    hist = pd.read_csv(csv_path)
    hist["table_type"] = hist["table_type"].astype(str).str.strip().str.lower()
    hist["phase"] = hist["phase"].astype(str).str.strip()
    hist["load_factor"] = hist["load_factor"].astype(float)
    return hist


def plot_probe_quantiles(hist: pd.DataFrame) -> None:
    if hist is None or hist.empty:
        return
    # one row per bin; p50/p90/p99 are constant per (cap, table, load, phase, repetition)
    qdf = (
        hist.groupby(["capacity_log2", "table_type", "load_factor", "phase", "repetition"])
        .agg(p50=("p50", "first"), p90=("p90", "first"), p99=("p99", "first"))
        .reset_index()
    )
    agg_q = (
        qdf.groupby(["capacity_log2", "table_type", "load_factor", "phase"])
        .agg(
            p50_mean=("p50", "mean"), p50_std=("p50", "std"),
            p90_mean=("p90", "mean"), p90_std=("p90", "std"),
            p99_mean=("p99", "mean"), p99_std=("p99", "std"),
        )
        .reset_index()
    )
    phases = agg_q["phase"].unique()
    capacities = sorted(agg_q["capacity_log2"].unique())
    for phase in phases:
        phase_df = agg_q[agg_q["phase"] == phase]
        fig, ax = plt.subplots(figsize=(7, 4))
        for cap in capacities:
            cap_df = phase_df[phase_df["capacity_log2"] == cap]
            for tbl in TABLE_ORDER:
                sub = cap_df[cap_df["table_type"] == tbl].sort_values("load_factor")
                if sub.empty:
                    continue
                for q, style in [("p50", "-"), ("p90", "--"), ("p99", ":")]:
                    ax.plot(
                        sub["load_factor"],
                        sub[f"{q}_mean"],
                        label=f"{tbl} (2^{cap}) {q}",
                        color=PALETTE.get(tbl, None),
                        linestyle=style,
                    )
        ax.set_xlabel("Load factor")
        ax.set_ylabel("Probe length (quantiles)")
        ax.set_title(f"Probe quantiles vs Load Factor: {phase}")
        ax.legend(title="Table / capacity / quantile", ncol=2, fontsize="small")
        fig.tight_layout()
        fig.savefig(f"probe_quantiles_{phase}.png")
        plt.close(fig)


def main() -> None:
    results_path = Path("results.csv")
    if not results_path.exists():
        raise SystemExit("results.csv not found in current directory.")

    df = load_results(results_path)
    agg = aggregate(df)

    plot_metric_vs_load(
        agg,
        metric_mean="time_per_op_mean",
        metric_std="time_per_op_std",
        ylabel="Time per operation (s)",
        prefix="time_per_op_vs_load",
        title_prefix="Time per operation vs Load Factor",
    )
    plot_metric_vs_load(
        agg,
        metric_mean="avg_probes_mean",
        metric_std="avg_probes_std",
        ylabel="Average probes per operation",
        prefix="probes_vs_load",
        title_prefix="Average Probes vs Load Factor",
    )
    plot_metric_vs_load(
        agg,
        metric_mean="max_probes_mean",
        metric_std="max_probes_std",
        ylabel="Max probes per operation",
        prefix="max_probes_vs_load",
        title_prefix="Max Probes vs Load Factor",
    )

    plot_time_vs_probes(agg)
    plot_boxplot_lookup_fail(df)

    hist_df = load_probe_hist(Path("probe_hist.csv"))
    plot_probe_quantiles(hist_df)

    # Console summary
    summary_cols = [
        "capacity_log2", "table_type", "load_factor", "phase",
        "time_per_op_mean", "time_per_op_std",
        "avg_probes_mean", "avg_probes_std",
        "max_probes_mean", "max_probes_std",
        "reps",
    ]
    print("\nAggregated metrics (means and std devs):")
    print(agg[summary_cols].to_string(index=False))
    print("\nPlots saved to current directory.")


if __name__ == "__main__":
    main()
