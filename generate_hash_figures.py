import os
import re
import math
import glob
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt


# ====================== CONFIG ======================

# Folder where all lp_*.csv and rh_*.csv live.
DATA_DIR = "output"      # change to "." if CSVs are in current dir

# Folder where figures will be saved.
FIG_DIR = "figures"

# Folder for summary CSVs.
SUMMARY_DIR = "summaries"


# =================== UTILITY STUFF ===================

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


# Example filenames:
#   lp_lf_25_insert_dynamic.csv
#   rh_lf_95_lookup_hit_dynamic.csv
#   rh_lf_95_final_distance_dynamic.csv
FILENAME_RE = re.compile(
    r'^(lp|rh)_lf_(\d+)_([a-z_]+)_(dynamic|fixed)\.csv$'
)
# groups:
#   1 -> scheme ('lp' or 'rh')
#   2 -> load factor as integer percent (e.g., "95")
#   3 -> metric type ("insert", "lookup", "lookup_hit", "lookup_miss", "final_distance", ...)
#   4 -> hist type ("dynamic" or "fixed")


def parse_filename(filename: str):
    """
    Parse a filename like 'lp_lf_95_lookup_hit_dynamic.csv'
      -> ('lp', 0.95, 'lookup_hit', 'dynamic')
    """
    base = os.path.basename(filename)
    m = FILENAME_RE.match(base)
    if not m:
        return None
    scheme = m.group(1)              # 'lp' or 'rh'
    lf_percent = int(m.group(2))     # 95
    metric = m.group(3)              # 'lookup_hit', 'lookup_miss', 'insert', ...
    hist_type = m.group(4)           # 'dynamic' or 'fixed'
    lf = lf_percent / 100.0
    return scheme, lf, metric, hist_type


def load_histograms():
    """
    Find all *_dynamic.csv and *_fixed.csv files, parse them, and organize as:

    dynamic_data[(scheme, metric)][load_factor] = DataFrame(probes,count)
      for all *dynamic.csv

    fixed_data[(scheme, metric)][load_factor]   = DataFrame(probes,count)
      for all *fixed.csv

    We mainly use dynamic_data for analysis, fixed_data is available if needed.
    """
    dynamic_data = defaultdict(dict)
    fixed_data = defaultdict(dict)

    pattern = os.path.join(DATA_DIR, "*.csv")
    for path in glob.glob(pattern):
        parsed = parse_filename(path)
        if parsed is None:
            continue
        scheme, lf, metric, hist_type = parsed
        df = pd.read_csv(path)
        # make sure it's sorted by probes
        df = df.sort_values("probes").reset_index(drop=True)

        if hist_type == "dynamic":
            dynamic_data[(scheme, metric)][lf] = df
        else:
            fixed_data[(scheme, metric)][lf] = df

    return dynamic_data, fixed_data


def hist_stats(df: pd.DataFrame):
    """
    Given a histogram DataFrame with columns: probes, count
    compute mean, variance, max, p90, p99.
    """
    probes = df["probes"].to_numpy()
    counts = df["count"].to_numpy()
    total = counts.sum()
    if total == 0:
        return dict(mean=0.0, var=0.0, max_probe=0, p90=0, p99=0)

    mean = (probes * counts).sum() / total
    mean_sq = ((probes ** 2) * counts).sum() / total
    var = max(0.0, mean_sq - mean ** 2)
    max_probe = int(probes.max())

    def percentile(p: float):
        """Discrete percentile based on histogram (ceil)."""
        target = math.ceil(p * total)
        cum = 0
        for k, c in zip(probes, counts):
            cum += c
            if cum >= target:
                return int(k)
        return int(probes.max())

    p90 = percentile(0.90)
    p99 = percentile(0.99)

    return dict(mean=mean, var=var, max_probe=max_probe, p90=p90, p99=p99)


def build_summary_table(dynamic_data, metric: str):
    """
    Build a summary DataFrame for a given metric, e.g. 'lookup_hit':
      columns: scheme, load_factor, mean, var, max_probe, p90, p99
    """
    rows = []
    for (scheme, m), by_lf in dynamic_data.items():
        if m != metric:
            continue
        for lf, df in by_lf.items():
            stats = hist_stats(df)
            row = {
                "scheme": scheme,
                "load_factor": lf,
                "metric": metric,
            }
            row.update(stats)
            rows.append(row)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.sort_values(["scheme", "load_factor"]).reset_index(drop=True)
    return df


def ccdf_from_hist(df: pd.DataFrame):
    """
    Compute CCDF: P(X >= x) from histogram DataFrame.
    Returns (probes, ccdf arrays), sorted by probes.
    """
    df_sorted = df.sort_values("probes")
    probes = df_sorted["probes"].to_numpy()
    counts = df_sorted["count"].to_numpy()
    total = counts.sum()
    if total == 0:
        return probes, counts

    # cumulative from right to left
    tail_counts = counts[::-1].cumsum()[::-1]
    ccdf = tail_counts / total
    return probes, ccdf


# =================== PLOTTING HELPERS ===================

def plot_stat_vs_load(summary_lp, summary_rh, stat_key: str, metric_label: str,
                      filename: str, logy=False):
    """
    Plot a stat (mean, p99, max_probe) vs load_factor
    for LP and RH summary tables for the same metric.
    """
    plt.figure()
    any_data = False

    if summary_lp is not None and not summary_lp.empty:
        plt.plot(summary_lp["load_factor"], summary_lp[stat_key],
                 marker="o", label="Linear Probing")
        any_data = True

    if summary_rh is not None and not summary_rh.empty:
        plt.plot(summary_rh["load_factor"], summary_rh[stat_key],
                 marker="s", label="Robin Hood")
        any_data = True

    if not any_data:
        plt.close()
        return

    plt.xlabel("Load factor")
    plt.ylabel(f"{metric_label} {stat_key}")
    plt.title(f"{metric_label}: {stat_key} vs load factor")
    if logy:
        plt.yscale("log")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    ensure_dir(FIG_DIR)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, filename))
    plt.close()


def plot_distribution_and_ccdf(dynamic_data, metric: str, lf: float,
                               filename_prefix: str, title_suffix: str = ""):
    """
    For a given metric (e.g. 'lookup_miss') and load factor lf, plot:
      - PMF (probe distribution) LP vs RH
      - CCDF LP vs RH
    """
    # PMF
    plt.figure()
    for scheme, label in [("lp", "Linear Probing"), ("rh", "Robin Hood")]:
        df = dynamic_data.get((scheme, metric), {}).get(lf)
        if df is None or df.empty:
            continue
        df_sorted = df.sort_values("probes")
        probes = df_sorted["probes"]
        probs = df_sorted["count"] / df_sorted["count"].sum()
        plt.step(probes, probs, where="mid", label=label)

    if not plt.gca().has_data():
        plt.close()
    else:
        plt.xlabel("Probes")
        plt.ylabel("Probability")
        title = f"{metric} probe distribution at load factor {lf:.2f}"
        if title_suffix:
            title += f" ({title_suffix})"
        plt.title(title)
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.legend()
        ensure_dir(FIG_DIR)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, f"{filename_prefix}_pmf_lf_{int(lf*100)}.png"))
        plt.close()

    # CCDF
    plt.figure()
    for scheme, label in [("lp", "Linear Probing"), ("rh", "Robin Hood")]:
        df = dynamic_data.get((scheme, metric), {}).get(lf)
        if df is None or df.empty:
            continue
        probes, ccdf = ccdf_from_hist(df)
        plt.step(probes, ccdf, where="post", label=label)

    if not plt.gca().has_data():
        plt.close()
    else:
        plt.xlabel("Probes")
        plt.ylabel("P(X ≥ x)")
        plt.yscale("log")
        title = f"{metric} CCDF at load factor {lf:.2f}"
        if title_suffix:
            title += f" ({title_suffix})"
        plt.title(title)
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, f"{filename_prefix}_ccdf_lf_{int(lf*100)}.png"))
        plt.close()


def plot_ccdf_ratio(dynamic_data, metric: str, lf: float, filename_prefix: str):
    """
    Plot ratio of CCDFs: P_LP(X >= x) / P_RH(X >= x) at a given lf.
    """
    df_lp = dynamic_data.get(("lp", metric), {}).get(lf)
    df_rh = dynamic_data.get(("rh", metric), {}).get(lf)
    if df_lp is None or df_lp.empty or df_rh is None or df_rh.empty:
        return

    probes_lp, ccdf_lp = ccdf_from_hist(df_lp)
    probes_rh, ccdf_rh = ccdf_from_hist(df_rh)

    # Align on common probe grid (union of probes)
    all_probes = sorted(set(probes_lp) | set(probes_rh))
    ccdf_lp_interp = []
    ccdf_rh_interp = []

    # Helper to get CCDF value for any k (piecewise constant)
    def ccdf_at(probes_arr, ccdf_arr, k):
        # probes_arr is sorted ascending
        # find index of first probe >= k
        # then return ccdf at that index or 0 if beyond range
        import bisect
        i = bisect.bisect_left(probes_arr, k)
        if i >= len(probes_arr):
            return 0.0
        return float(ccdf_arr[i])

    for k in all_probes:
        lp_val = ccdf_at(probes_lp, ccdf_lp, k)
        rh_val = ccdf_at(probes_rh, ccdf_rh, k)
        ccdf_lp_interp.append(lp_val)
        ccdf_rh_interp.append(rh_val)

    ratio = []
    for lp, rh in zip(ccdf_lp_interp, ccdf_rh_interp):
        if rh == 0:
            ratio.append(float("inf") if lp > 0 else 1.0)
        else:
            ratio.append(lp / rh)

    # Now plot ratio
    plt.figure()
    plt.step(all_probes, ratio, where="post")
    plt.xlabel("Probes")
    plt.ylabel("CCDF ratio LP/RH")
    plt.yscale("log")
    plt.title(f"{metric}: CCDF ratio LP/RH at load factor {lf:.2f}")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    ensure_dir(FIG_DIR)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f"{filename_prefix}_ccdf_ratio_lf_{int(lf*100)}.png"))
    plt.close()


# ========================= MAIN =========================

def main():
    ensure_dir(FIG_DIR)
    ensure_dir(SUMMARY_DIR)

    dynamic_data, fixed_data = load_histograms()

    # Metrics we know about / care about
    metrics_of_interest = [
        "insert",
        "lookup",
        "lookup_hit",
        "lookup_miss",
        "final_distance",   # RH-only
    ]

    # Build summary tables for each metric
    summaries = {}
    for metric in metrics_of_interest:
        summaries[metric] = build_summary_table(dynamic_data, metric=metric)

    # Write combined summary to CSV
    all_rows = []
    for metric, df in summaries.items():
        if df.empty:
            continue
        all_rows.append(df)
    if all_rows:
        summary_all = pd.concat(all_rows, ignore_index=True)
        summary_all = summary_all.sort_values(["metric", "scheme", "load_factor"])
        summary_all.to_csv(os.path.join(SUMMARY_DIR, "summary_all_metrics.csv"), index=False)

    # Convenience: split LP vs RH for each metric
    def split_scheme(metric_name):
        df = summaries.get(metric_name)
        if df is None or df.empty:
            return None, None
        return (df[df["scheme"] == "lp"].copy(),
                df[df["scheme"] == "rh"].copy())

    # ---- 1) Stats vs load factor for each metric ----

    # Insert (LP vs RH)
    insert_lp, insert_rh = split_scheme("insert")
    plot_stat_vs_load(insert_lp, insert_rh, "mean", "Insert", "insert_mean_vs_load.png", logy=False)
    plot_stat_vs_load(insert_lp, insert_rh, "p99", "Insert", "insert_p99_vs_load.png", logy=True)
    plot_stat_vs_load(insert_lp, insert_rh, "max_probe", "Insert", "insert_max_vs_load.png", logy=True)

    # Lookup (aggregated)
    lookup_lp, lookup_rh = split_scheme("lookup")
    plot_stat_vs_load(lookup_lp, lookup_rh, "mean", "Lookup (all)", "lookup_mean_vs_load.png", logy=False)
    plot_stat_vs_load(lookup_lp, lookup_rh, "p99", "Lookup (all)", "lookup_p99_vs_load.png", logy=True)
    plot_stat_vs_load(lookup_lp, lookup_rh, "max_probe", "Lookup (all)", "lookup_max_vs_load.png", logy=True)

    # Lookup HIT
    hit_lp, hit_rh = split_scheme("lookup_hit")
    plot_stat_vs_load(hit_lp, hit_rh, "mean", "Lookup HIT", "lookup_hit_mean_vs_load.png", logy=False)
    plot_stat_vs_load(hit_lp, hit_rh, "p99", "Lookup HIT", "lookup_hit_p99_vs_load.png", logy=True)
    plot_stat_vs_load(hit_lp, hit_rh, "max_probe", "Lookup HIT", "lookup_hit_max_vs_load.png", logy=True)

    # Lookup MISS
    miss_lp, miss_rh = split_scheme("lookup_miss")
    plot_stat_vs_load(miss_lp, miss_rh, "mean", "Lookup MISS", "lookup_miss_mean_vs_load.png", logy=False)
    plot_stat_vs_load(miss_lp, miss_rh, "p99", "Lookup MISS", "lookup_miss_p99_vs_load.png", logy=True)
    plot_stat_vs_load(miss_lp, miss_rh, "max_probe", "Lookup MISS", "lookup_miss_max_vs_load.png", logy=True)

    # Final distance (RH only)
    final_lp, final_rh = split_scheme("final_distance")
    # LP will be empty, RH-only
    plot_stat_vs_load(None, final_rh, "mean", "RH final displacement", "rh_final_distance_mean_vs_load.png", logy=False)
    plot_stat_vs_load(None, final_rh, "p99", "RH final displacement", "rh_final_distance_p99_vs_load.png", logy=True)
    plot_stat_vs_load(None, final_rh, "max_probe", "RH final displacement", "rh_final_distance_max_vs_load.png", logy=True)

    # ---- 2) Distribution + CCDF at highest load factor ----

    # Choose highest LF where we have lookup_miss LP entries
    highest_lf = None
    if miss_lp is not None and not miss_lp.empty:
        highest_lf = miss_lp["load_factor"].max()
    elif lookup_lp is not None and not lookup_lp.empty:
        highest_lf = lookup_lp["load_factor"].max()

    if highest_lf is not None:
        # At highest LF: distributions and CCDFs for lookup_hit, lookup_miss, lookup(all)
        plot_distribution_and_ccdf(dynamic_data, "lookup_hit", highest_lf,
                                   filename_prefix="lookup_hit",
                                   title_suffix="highest load factor")
        plot_distribution_and_ccdf(dynamic_data, "lookup_miss", highest_lf,
                                   filename_prefix="lookup_miss",
                                   title_suffix="highest load factor")
        plot_distribution_and_ccdf(dynamic_data, "lookup", highest_lf,
                                   filename_prefix="lookup_all",
                                   title_suffix="highest load factor")

        # CCDF ratio LP/RH for lookup_miss at highest LF (most interesting tail)
        plot_ccdf_ratio(dynamic_data, "lookup_miss", highest_lf,
                        filename_prefix="lookup_miss")

        # RH final displacement distribution at highest LF (if present)
        rh_finals = dynamic_data.get(("rh", "final_distance"), {})
        if rh_finals:
            # pick max LF where we have final_distance
            lf_fd = max(rh_finals.keys())
            plot_distribution_and_ccdf(dynamic_data, "final_distance", lf_fd,
                                       filename_prefix="rh_final_distance",
                                       title_suffix="final layout (Robin Hood only)")

    print("Done. Figures in:", FIG_DIR)
    print("Summary CSVs in:", SUMMARY_DIR)


if __name__ == "__main__":
    main()
