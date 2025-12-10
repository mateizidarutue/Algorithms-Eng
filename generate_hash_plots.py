import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

# ================= CONFIGURATION =================
# Set to the folder containing the CSV files
ROOT_DIR = 'output' 
OUTPUT_DIR = 'research_plots_v2'  # new folder to avoid mixing with old plots
# =================================================

def parse_filename_flexible(filename):
    """
    Parses filenames to extract Algorithm, Metric, Type, and Load Factor.
    Expected components in filename:
      - 'lp' or 'rh' (Linear Probing vs Robin Hood)
      - 'insert' or 'lookup' or 'final_distance'
      - 'fixed' or 'dynamic'
      - 'lf_XX' (e.g., lf_0.99, lf_0.5)
    """
    name = filename.lower()
    
    # 1. Detect Algorithm
    if 'lp' in name and 'rh' not in name:
        algo = 'lp'
    elif 'rh' in name and 'lp' not in name:
        algo = 'rh'
    else:
        return None # Can't determine algo or ambiguous

    # 2. Detect Metric
    if 'lookup' in name:
        metric = 'lookup'
    elif 'insert' in name:
        metric = 'insert'
    elif 'final_distance' in name:
        metric = 'final_distance'
    else:
        return None

    # 3. Detect Type (Fixed vs Dynamic)
    if 'dynamic' in name:
        ftype = 'dynamic'
    elif 'fixed' in name:
        ftype = 'fixed'
    else:
        return None

    # 4. Extract Load Factor (lf_0.99 or lf_99)
    # Looks for "lf_" followed by numbers/dots
    lf_match = re.search(r'lf_([0-9\.]+)', name)
    if lf_match:
        try:
            lf_val = float(lf_match.group(1))
            lf = lf_val if lf_val <= 1.5 else lf_val / 100.0
        except ValueError:
            return None
    else:
        return None

    return {
        'algo': algo,
        'metric': metric,
        'type': ftype,
        'lf': lf
    }

def load_csv(filepath):
    if not filepath.exists(): return None
    try:
        df = pd.read_csv(filepath)
        if 'probes' in df.columns and 'count' in df.columns:
            return df.sort_values('probes')
    except Exception:
        pass
    return None

def add_file(groups, meta, path):
    lf = meta['lf']
    algo = meta['algo']
    metric = meta['metric']
    ftype = meta['type']
    groups.setdefault(lf, {}).setdefault(algo, {}).setdefault(metric, {})[ftype] = path

def select_histogram(groups, lf, algo, metric, preferred='dynamic'):
    algo_group = groups.get(lf, {}).get(algo, {}).get(metric, {})
    return algo_group.get(preferred) or algo_group.get('dynamic') or algo_group.get('fixed')

def prepare_histogram(df, max_bins=200):
    """Reduce extremely wide histograms by binning when there are too many distinct probe counts."""
    if df is None or df.empty:
        return df, 1
    unique_bins = df['probes'].nunique()
    if unique_bins <= max_bins:
        return df, 1
    max_probe = int(df['probes'].max())
    bin_width = max(1, int(np.ceil(max_probe / max_bins)))
    tmp = df.copy()
    tmp['probe_bin'] = (tmp['probes'] // bin_width) * bin_width
    collapsed = tmp.groupby('probe_bin')['count'].sum().reset_index().rename(columns={'probe_bin': 'probes'})
    collapsed = collapsed.sort_values('probes')
    return collapsed, bin_width

def build_distribution(df):
    """Return DataFrame with probes, count, prob, cdf, ccdf sorted by probes."""
    if df is None or df.empty:
        return None
    df = df.sort_values('probes').copy()
    total = df['count'].sum()
    if total == 0:
        return None
    df['prob'] = df['count'] / total
    df['cdf'] = df['prob'].cumsum()
    df['ccdf'] = 1.0 - df['cdf'] + df['prob']  # include the current bin
    return df[['probes', 'count', 'prob', 'cdf', 'ccdf']]

def calculate_weighted_stats(df, tail_thresholds=(2, 4, 8, 16)):
    if df is None or df.empty: return {}
    
    total = df['count'].sum()
    if total == 0: return {}

    weighted_sum = (df['probes'] * df['count']).sum()
    avg = weighted_sum / total
    
    weighted_sq_sum = ((df['probes'] ** 2) * df['count']).sum()
    variance = (weighted_sq_sum / total) - (avg ** 2)
    
    cumsum = df['count'].cumsum()
    def get_percentile(p):
        idx = cumsum.searchsorted(p * total)
        return float(df.iloc[min(idx, len(df)-1)]['probes'])

    stats = {
        'avg': avg,
        'var': variance,
        'max': float(df['probes'].max()),
        'p50': get_percentile(0.5),
        'p90': get_percentile(0.9),
        'p99': get_percentile(0.99),
    }

    for t in tail_thresholds:
        stats[f'ge_{t}'] = float(df[df['probes'] >= t]['count'].sum() / total)

    return stats

def plot_comparative_cdf(lp_df, rh_df, load_factor, metric, output_path):
    plt.figure(figsize=(10, 6))
    datasets = [('Linear Probing', lp_df, 'blue'), ('Robin Hood', rh_df, 'orange')]
    has_zero = False
    for name, df, color in datasets:
        if df is None: continue
        df = df.copy()
        has_zero = has_zero or (df['probes'] == 0).any()
        total = df['count'].sum()
        df['cdf'] = df['count'].cumsum() / total
        plt.plot(df['probes'], df['cdf'], label=name, color=color, linewidth=2.5)

    plt.title(f"CDF of {metric.title()} Probes (Load Factor {load_factor:.2f})", fontsize=14)
    plt.xlabel("Probe Count", fontsize=12)
    plt.ylabel("Probability", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(fontsize=12)
    if not has_zero:
        plt.xscale('log')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_ccdf(lp_df, rh_df, lp_stats, rh_stats, load_factor, metric, output_path):
    """Complementary CDF with log-y and quantile markers to highlight tails clearly."""
    lp_dist = build_distribution(lp_df)
    rh_dist = build_distribution(rh_df)
    if lp_dist is None or rh_dist is None:
        return

    plt.figure(figsize=(10, 6))
    for name, dist, color, stats in [
        ('Linear Probing', lp_dist, 'blue', lp_stats),
        ('Robin Hood', rh_dist, 'orange', rh_stats)
    ]:
        plt.step(dist['probes'], dist['ccdf'], where='post', label=name, color=color, linewidth=2.2)
        for q in ['p50', 'p90', 'p99']:
            if q in stats and stats[q] is not None:
                plt.axvline(stats[q], color=color, linestyle='--', alpha=0.35)

    plt.yscale('log')
    plt.xlabel('Probes', fontsize=12)
    plt.ylabel('P(Probes ≥ x)', fontsize=12)
    plt.title(f'CCDF of {metric.title()} Probes (LF={load_factor:.2f})', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_pmf_step(lp_df, rh_df, load_factor, metric, output_path):
    """Probability mass (normalized) with log-y to make both distributions visible."""
    lp_dist = build_distribution(lp_df)
    rh_dist = build_distribution(rh_df)
    if lp_dist is None or rh_dist is None:
        return
    plt.figure(figsize=(10, 6))
    plt.step(lp_dist['probes'], lp_dist['prob'], where='mid', label='Linear Probing', color='blue', linewidth=2.0)
    plt.step(rh_dist['probes'], rh_dist['prob'], where='mid', label='Robin Hood', color='orange', linewidth=2.0)
    plt.yscale('log')
    plt.xlabel('Probes', fontsize=12)
    plt.ylabel('Probability mass', fontsize=12)
    plt.title(f'Probe PMF (normalized) - {metric.title()} (LF={load_factor:.2f})', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_ccdf_ratio(lp_df, rh_df, load_factor, metric, output_path, epsilon=1e-12):
    """Ratio of CCDFs (LP / RH) to show multiplicative tail advantage."""
    lp_dist = build_distribution(lp_df)
    rh_dist = build_distribution(rh_df)
    if lp_dist is None or rh_dist is None:
        return

    x_vals = sorted(set(lp_dist['probes']).union(set(rh_dist['probes'])))
    x_arr = np.array(x_vals)

    def values_at(dist):
        probes = dist['probes'].to_numpy()
        ccdf_vals = dist['ccdf'].to_numpy()
        idx = np.searchsorted(probes, x_arr, side='right') - 1
        idx = np.clip(idx, -1, len(probes) - 1)
        vals = np.where(idx >= 0, ccdf_vals[idx], 1.0)
        return vals

    lp_ccdf = values_at(lp_dist)
    rh_ccdf = values_at(rh_dist)
    ratio = (lp_ccdf + epsilon) / (rh_ccdf + epsilon)

    plt.figure(figsize=(10, 5))
    plt.plot(x_arr, ratio, color='purple', linewidth=2.0)
    plt.axhline(1.0, color='k', linestyle='--', alpha=0.6)
    plt.xlabel('Probes', fontsize=12)
    plt.ylabel('LP CCDF / RH CCDF', fontsize=12)
    plt.title(f'Tail Advantage (CCDF ratio) - {metric.title()} (LF={load_factor:.2f})', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_trends(agg_df, metric, output_dir):
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return

    sns.set_style("whitegrid")
    
    metrics = [
        ('lp_max', 'rh_max', 'Worst-Case Probes', 'Max Probes', f'{metric}_trend_worst_case.png'),
        ('lp_p99', 'rh_p99', 'Tail Latency (P99)', '99th % Probes', f'{metric}_trend_p99.png'),
        ('lp_var', 'rh_var', 'Probe Variance', 'Variance', f'{metric}_trend_variance.png')
    ]
    
    for lp_col, rh_col, title, ylabel, fname in metrics:
        plt.figure(figsize=(8, 5))
        plt.plot(df['load_factor'], df[lp_col], 'o--', label='Linear Probing', color='blue')
        plt.plot(df['load_factor'], df[rh_col], 's-', label='Robin Hood', color='orange')
        plt.xlabel('Load Factor')
        plt.ylabel(ylabel)
        plt.title(f"{title} ({metric.title()})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir / fname, dpi=300)
        plt.close()

def plot_histogram_comparison(lp_df, rh_df, load_factor, metric, output_path):
    """Compare probe count distributions as histograms with side-by-side bars"""
    fig, ax = plt.subplots(figsize=(14, 6))

    lp_df, lp_bin = prepare_histogram(lp_df)
    rh_df, rh_bin = prepare_histogram(rh_df)
    bin_hint = f" (bin width {max(lp_bin, rh_bin)})" if max(lp_bin, rh_bin) > 1 else ""
    
    if lp_df is not None and not lp_df.empty and rh_df is not None and not rh_df.empty:
        # Get probe counts from both dataframes
        lp_probes = set(lp_df['probes'].values)
        rh_probes = set(rh_df['probes'].values)
        all_probes = sorted(lp_probes.union(rh_probes))
        
        # Create dictionaries for easy lookup
        lp_counts = dict(zip(lp_df['probes'], lp_df['count']))
        rh_counts = dict(zip(rh_df['probes'], rh_df['count']))
        
        # Prepare data
        lp_vals = [lp_counts.get(p, 0) for p in all_probes]
        rh_vals = [rh_counts.get(p, 0) for p in all_probes]
        
        # Create side-by-side bars
        x = range(len(all_probes))
        width = 0.35
        
        bars1 = ax.bar([i - width/2 for i in x], lp_vals, width, label='Linear Probing', color='#1f77b4', alpha=0.8)
        bars2 = ax.bar([i + width/2 for i in x], rh_vals, width, label='Robin Hood', color='#ff7f0e', alpha=0.8)
        
        ax.set_xlabel('Probe Count', fontsize=12, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax.set_title(f"Probe Distribution: {metric.upper()} (Load Factor {load_factor:.2f}){bin_hint}", fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(all_probes, rotation=45 if len(all_probes) > 20 else 0)
        ax.legend(fontsize=11, loc='upper right')
        ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    else:
        # Fallback if one dataframe is missing
        if lp_df is not None and not lp_df.empty:
            ax.bar(lp_df['probes'], lp_df['count'], alpha=0.8, label='Linear Probing', color='#1f77b4', width=0.6)
        if rh_df is not None and not rh_df.empty:
            ax.bar(rh_df['probes'], rh_df['count'], alpha=0.8, label='Robin Hood', color='#ff7f0e', width=0.6)
        ax.set_xlabel('Probe Count', fontsize=12, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax.set_title(f"Probe Distribution: {metric.upper()} (Load Factor {load_factor:.2f}){bin_hint}", fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_avg_comparison(agg_df, metric, output_dir):
    """Compare average probe counts across load factors"""
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return
    
    plt.figure(figsize=(10, 6))
    plt.plot(df['load_factor'], df['lp_avg'], 'o--', label='Linear Probing', color='blue', linewidth=2.5, markersize=8)
    plt.plot(df['load_factor'], df['rh_avg'], 's-', label='Robin Hood', color='orange', linewidth=2.5, markersize=8)
    plt.xlabel('Load Factor', fontsize=12)
    plt.ylabel('Average Probes', fontsize=12)
    plt.title(f'Average Probe Count vs Load Factor ({metric.title()})', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_trend_avg_probes.png', dpi=300)
    plt.close()

def plot_ratio_improvement(agg_df, metric, output_dir):
    """Show how much better Robin Hood is compared to Linear Probing"""
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return
    
    # Calculate improvement ratios
    df['var_ratio'] = df['lp_var'] / (df['rh_var'] + 1e-10)  # avoid division by zero
    df['max_ratio'] = df['lp_max'] / (df['rh_max'] + 1e-10)
    df['p99_ratio'] = df['lp_p99'] / (df['rh_p99'] + 1e-10)
    df['avg_ratio'] = df['lp_avg'] / (df['rh_avg'] + 1e-10)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    plot_meta = [
        ('max_ratio', 'Max Probes'),
        ('p99_ratio', 'P99 Probes'),
        ('var_ratio', 'Variance'),
        ('avg_ratio', 'Average Probes')
    ]
    
    for ax, (col, title) in zip(axes.flat, plot_meta):
        ax.plot(df['load_factor'], df[col], 'o-', linewidth=2.0)
        ax.set_title(f'{title} (LP / RH)', fontsize=11)
        ax.set_xlabel('Load Factor')
        ax.set_ylabel('Improvement Ratio')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.5)
    
    plt.suptitle(f'Improvement Ratios ({metric.title()})', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_dir / f'{metric}_improvement_ratios.png', dpi=300)
    plt.close()

def plot_quantile_panels(agg_df, metric, output_dir):
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for ax, q in zip(axes, ['p50', 'p90', 'p99']):
        ax.plot(df['load_factor'], df[f'lp_{q}'], 'o--', label='Linear Probing', color='blue')
        ax.plot(df['load_factor'], df[f'rh_{q}'], 's-', label='Robin Hood', color='orange')
        ax.set_xlabel('Load Factor')
        ax.set_ylabel(f'{q.upper()} Probes')
        ax.set_title(f'{q.upper()} vs Load ({metric.title()})')
        ax.grid(True, linestyle='--', alpha=0.3)
        if q == 'p50':
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_trend_quantiles.png', dpi=300)
    plt.close()

def plot_tail_probabilities(agg_df, metric, thresholds, output_dir):
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return

    fig, axes = plt.subplots(1, len(thresholds), figsize=(5 * len(thresholds), 4), sharey=True)
    for ax, t in zip(axes, thresholds):
        ax.plot(df['load_factor'], df[f'lp_ge_{t}'], 'o--', label='Linear Probing', color='blue')
        ax.plot(df['load_factor'], df[f'rh_ge_{t}'], 's-', label='Robin Hood', color='orange')
        ax.set_xlabel('Load Factor')
        ax.set_ylabel(f'P(probes >= {t})')
        ax.set_title(f'Tail Probability >= {t} ({metric.title()})')
        ax.grid(True, linestyle='--', alpha=0.3)
        if t == thresholds[0]:
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_tail_probabilities.png', dpi=300)
    plt.close()

def plot_summary_absolute(agg_df, metric, output_dir):
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return
    stats = [('avg', 'Average'), ('p90', 'P90'), ('p99', 'P99'), ('ge_8', 'P≥8'), ('ge_16', 'P≥16')]
    cols = 3
    rows = int(np.ceil(len(stats) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.2, rows * 3.4))
    axes = axes.flatten()
    for ax, (key, label) in zip(axes, stats):
        ax.plot(df['load_factor'], df[f'lp_{key}'], 'o--', label='LP', color='blue')
        ax.plot(df['load_factor'], df[f'rh_{key}'], 's-', label='RH', color='orange')
        ax.set_xlabel('Load Factor')
        ax.set_ylabel(label)
        ax.set_title(label)
        ax.grid(True, alpha=0.3)
    axes[0].legend()
    plt.suptitle(f'{metric.title()} - Absolute Metrics', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_dir / f'{metric}_summary_absolute.png', dpi=300)
    plt.close()

def plot_summary_ratio(agg_df, metric, output_dir, epsilon=1e-12):
    df = agg_df[agg_df['metric'] == metric].sort_values('load_factor')
    if df.empty: return
    stats = [('avg', 'Average'), ('p90', 'P90'), ('p99', 'P99'), ('ge_8', 'P≥8'), ('ge_16', 'P≥16')]
    cols = 3
    rows = int(np.ceil(len(stats) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.2, rows * 3.4))
    axes = axes.flatten()
    for ax, (key, label) in zip(axes, stats):
        ratio = (df[f'lp_{key}'] + epsilon) / (df[f'rh_{key}'] + epsilon)
        ax.plot(df['load_factor'], ratio, 'o-', color='purple')
        ax.axhline(1.0, color='k', linestyle='--', alpha=0.6)
        ax.set_xlabel('Load Factor')
        ax.set_ylabel('LP / RH')
        ax.set_title(f'{label} ratio')
        ax.grid(True, alpha=0.3)
    plt.suptitle(f'{metric.title()} - LP vs RH Ratios', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_dir / f'{metric}_summary_ratio.png', dpi=300)
    plt.close()

def plot_insert_vs_lookup(groups, output_dir, preferred_type='dynamic'):
    """Compare insert vs lookup performance"""
    plt.figure(figsize=(12, 6))
    
    lfs = []
    lp_insert_avgs = []
    lp_lookup_avgs = []
    rh_insert_avgs = []
    rh_lookup_avgs = []
    
    for lf in sorted(groups.keys()):
        lfs.append(lf)
        
        lp_insert = load_csv(select_histogram(groups, lf, 'lp', 'insert', preferred_type)) if select_histogram(groups, lf, 'lp', 'insert', preferred_type) else None
        lp_lookup = load_csv(select_histogram(groups, lf, 'lp', 'lookup', preferred_type)) if select_histogram(groups, lf, 'lp', 'lookup', preferred_type) else None
        rh_insert = load_csv(select_histogram(groups, lf, 'rh', 'insert', preferred_type)) if select_histogram(groups, lf, 'rh', 'insert', preferred_type) else None
        rh_lookup = load_csv(select_histogram(groups, lf, 'rh', 'lookup', preferred_type)) if select_histogram(groups, lf, 'rh', 'lookup', preferred_type) else None
        
        lp_insert_stats = calculate_weighted_stats(lp_insert)
        lp_lookup_stats = calculate_weighted_stats(lp_lookup)
        rh_insert_stats = calculate_weighted_stats(rh_insert)
        rh_lookup_stats = calculate_weighted_stats(rh_lookup)
        
        lp_insert_avgs.append(lp_insert_stats.get('avg', 0))
        lp_lookup_avgs.append(lp_lookup_stats.get('avg', 0))
        rh_insert_avgs.append(rh_insert_stats.get('avg', 0))
        rh_lookup_avgs.append(rh_lookup_stats.get('avg', 0))

    if not lfs:
        return
    
    x = range(len(lfs))
    width = 0.2
    
    plt.bar([i - 1.5*width for i in x], lp_insert_avgs, width, label='LP Insert', color='blue', alpha=0.7)
    plt.bar([i - 0.5*width for i in x], lp_lookup_avgs, width, label='LP Lookup', color='lightblue', alpha=0.7)
    plt.bar([i + 0.5*width for i in x], rh_insert_avgs, width, label='RH Insert', color='orange', alpha=0.7)
    plt.bar([i + 1.5*width for i in x], rh_lookup_avgs, width, label='RH Lookup', color='lightyellow', alpha=0.7)
    
    plt.xlabel('Load Factor', fontsize=12)
    plt.ylabel('Average Probes', fontsize=12)
    plt.title('Insert vs Lookup Performance Comparison', fontsize=14)
    plt.xticks(x, [f'{lf:.2f}' for lf in lfs])
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(output_dir / 'insert_vs_lookup.png', dpi=300)
    plt.close()

def main():
    root = Path(ROOT_DIR)
    out = Path(OUTPUT_DIR)
    out.mkdir(exist_ok=True)

    groups = {}
    tail_thresholds = (2, 4, 8, 16)
    preferred_hist = 'dynamic'

    print(f"Scanning directory: {root.resolve()}")
    
    # SCAN FILES
    for entry in root.glob("*.csv"):
        meta = parse_filename_flexible(entry.name)
        if not meta:
            continue # Skip unrecognized files
        add_file(groups, meta, entry)

    # PROCESS GROUPS
    agg_records = []
    sorted_lfs = sorted(groups.keys())
    
    if not sorted_lfs:
        print("No valid files found. Ensure filenames contain 'lp'/'rh', 'dynamic', and 'lf_XX'.")
        return

    print(f"Found Load Factors: {sorted_lfs}")

    for lf in sorted_lfs:
        print(f"  -> Generating plots for LF {lf:.2f}...")
        for metric in ['lookup', 'insert']:
            lp_path = select_histogram(groups, lf, 'lp', metric, preferred_hist)
            rh_path = select_histogram(groups, lf, 'rh', metric, preferred_hist)

            lp_df = load_csv(lp_path) if lp_path else None
            rh_df = load_csv(rh_path) if rh_path else None

            if lp_df is None or rh_df is None:
                print(f"     - Skipping {metric} (missing {'LP' if lp_df is None else ''}{' & ' if lp_df is None and rh_df is None else ''}{'RH' if rh_df is None else ''})")
                continue

            lp_stats = calculate_weighted_stats(lp_df, tail_thresholds)
            rh_stats = calculate_weighted_stats(rh_df, tail_thresholds)

            # Per-load-factor visuals that highlight both distributions
            plot_ccdf(lp_df, rh_df, lp_stats, rh_stats, lf, metric, out / f"ccdf_{metric}_lf_{lf:.2f}.png")
            plot_pmf_step(lp_df, rh_df, lf, metric, out / f"pmf_{metric}_lf_{lf:.2f}.png")
            plot_ccdf_ratio(lp_df, rh_df, lf, metric, out / f"ccdf_ratio_{metric}_lf_{lf:.2f}.png")

            record = {'load_factor': lf, 'metric': metric}
            for prefix, stats in [('lp', lp_stats), ('rh', rh_stats)]:
                for key, val in stats.items():
                    record[f"{prefix}_{key}"] = val
            agg_records.append(record)

    # GENERATE TRENDS
    agg_df = pd.DataFrame(agg_records)
    if not agg_df.empty:
        print("Generating Aggregate Trend Plots...")
        for metric in ['lookup', 'insert']:
            plot_trends(agg_df, metric, out)
            plot_avg_comparison(agg_df, metric, out)
            plot_ratio_improvement(agg_df, metric, out)
            plot_quantile_panels(agg_df, metric, out)
            plot_tail_probabilities(agg_df, metric, tail_thresholds, out)
            plot_summary_absolute(agg_df, metric, out)
            plot_summary_ratio(agg_df, metric, out)
        plot_insert_vs_lookup(groups, out, preferred_type=preferred_hist)
        print(f"Done. Visualizations saved to '{OUTPUT_DIR}'")
    else:
        print("Not enough paired data (LP + RH) to generate trends.")

if __name__ == "__main__":
    main()
