#!/usr/bin/env python3
"""
generate_ccdf_highload.py

Build CCDF plots for high-load factors for:
- Insert probes (LP vs RH)
- Lookup HIT probes (LP vs RH)
- Lookup MISS probes (LP vs RH)
- Robin Hood displacement from home slot

Uses the *_dynamic.csv histograms you exported.
"""

import csv
import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

# High-load factors to plot (as percentages used in filenames)
HIGH_LOAD_LFS = [80, 90, 95, 96, 97, 98, 99]

# Where CSVs live and where to put figures
BASE_DIR = Path("output")        # change if needed
FIG_DIR = Path("figs_ccdf_highload")  # output directory for figures
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Filename patterns; {scheme} is "lp" or "rh", {lf} is e.g. 80
FILENAME_TEMPLATES: Dict[str, str] = {
    "insert": "{scheme}_lf_{lf}_insert_dynamic.csv",
    "lookup_hit": "{scheme}_lf_{lf}_lookup_hit_dynamic.csv",
    "lookup_miss": "{scheme}_lf_{lf}_lookup_miss_dynamic.csv",
    # Robin Hood only:
    "rh_displacement": "rh_lf_{lf}_final_distance_dynamic.csv",
}

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

def load_histogram(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load a two-column CSV (probes,count) and return (values, counts)
    sorted by probe length.
    """
    xs: List[int] = []
    cs: List[int] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if "probes" not in reader.fieldnames or "count" not in reader.fieldnames:
            raise ValueError(f"{path} does not have 'probes,count' columns")
        for row in reader:
            xs.append(int(row["probes"]))
            cs.append(int(row["count"]))

    xs_arr = np.array(xs, dtype=np.int64)
    cs_arr = np.array(cs, dtype=np.int64)

    # sort by probe length
    order = np.argsort(xs_arr)
    return xs_arr[order], cs_arr[order]


def ccdf_from_hist(xs: np.ndarray, cs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Given probe lengths xs and their counts cs, compute CCDF:
        P(X >= x_i).
    Returns (xs, ccdf_values) with same length as xs.
    """
    total = cs.sum()
    if total == 0:
        return xs, np.zeros_like(xs, dtype=float)
    # cumulative from the right
    tail_counts = np.cumsum(cs[::-1])[::-1]
    ccdf_vals = tail_counts / total
    return xs, ccdf_vals


def plot_ccdf_pair(
    lf: int,
    metric: str,
    lp_hist: Tuple[np.ndarray, np.ndarray],
    rh_hist: Tuple[np.ndarray, np.ndarray],
):
    """
    Plot CCDF for LP and RH on the same figure.
    """
    lp_x, lp_c = lp_hist
    rh_x, rh_c = rh_hist
    lp_xc, lp_ccdf = ccdf_from_hist(lp_x, lp_c)
    rh_xc, rh_ccdf = ccdf_from_hist(rh_x, rh_c)

    plt.figure(figsize=(7, 4.5))
    plt.step(lp_xc, lp_ccdf, where="post", label="Linear Probing")
    plt.step(rh_xc, rh_ccdf, where="post", label="Robin Hood")

    plt.yscale("log")
    plt.xlabel("Probes")
    title_metric = {
        "insert": "Insert",
        "lookup_hit": "Lookup HIT",
        "lookup_miss": "Lookup MISS",
    }[metric]
    plt.title(f"{title_metric} probes CCDF at load factor {lf/100:.2f}")
    plt.ylabel("P(probes ≥ x)")
    plt.grid(True, which="both", linestyle=":", linewidth=0.5)
    plt.legend()

    fname = FIG_DIR / f"ccdf_{metric}_lf_{lf}.png"
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"Saved {fname}")


def plot_ccdf_rh_displacement(lf: int, hist: Tuple[np.ndarray, np.ndarray]):
    """
    CCDF for Robin Hood displacement from home slot.
    """
    xs, cs = hist
    xs_c, ccdf = ccdf_from_hist(xs, cs)

    plt.figure(figsize=(7, 4.5))
    plt.step(xs_c, ccdf, where="post", label="Robin Hood displacement")

    plt.yscale("log")
    plt.xlabel("Displacement from home slot")
    plt.ylabel("P(displacement ≥ x)")
    plt.title(f"Robin Hood displacement CCDF at load factor {lf/100:.2f}")
    plt.grid(True, which="both", linestyle=":", linewidth=0.5)
    plt.legend()

    fname = FIG_DIR / f"ccdf_rh_displacement_lf_{lf}.png"
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"Saved {fname}")


# -------------------------------------------------------------------
# Main driver
# -------------------------------------------------------------------

def main():
    for lf in HIGH_LOAD_LFS:
        print(f"\n=== Load factor {lf/100:.2f} ===")

        # Insert, HIT, MISS (LP vs RH)
        for metric in ["insert", "lookup_hit", "lookup_miss"]:
            try:
                lp_file = BASE_DIR / FILENAME_TEMPLATES[metric].format(
                    scheme="lp", lf=lf
                )
                rh_file = BASE_DIR / FILENAME_TEMPLATES[metric].format(
                    scheme="rh", lf=lf
                )

                if not lp_file.exists():
                    print(f"  [WARN] Missing {lp_file}, skipping {metric} for LF={lf}")
                    continue
                if not rh_file.exists():
                    print(f"  [WARN] Missing {rh_file}, skipping {metric} for LF={lf}")
                    continue

                lp_hist = load_histogram(lp_file)
                rh_hist = load_histogram(rh_file)
                plot_ccdf_pair(lf, metric, lp_hist, rh_hist)
            except Exception as e:
                print(f"  [ERROR] {metric} LF={lf}: {e}")

        # Robin Hood displacement only
        try:
            disp_template = FILENAME_TEMPLATES["rh_displacement"]
            rh_disp_file = BASE_DIR / disp_template.format(lf=lf)
            if rh_disp_file.exists():
                rh_disp_hist = load_histogram(rh_disp_file)
                plot_ccdf_rh_displacement(lf, rh_disp_hist)
            else:
                print(f"  [WARN] No RH displacement file for LF={lf}: {rh_disp_file}")
        except Exception as e:
            print(f"  [ERROR] RH displacement LF={lf}: {e}")


if __name__ == "__main__":
    main()
