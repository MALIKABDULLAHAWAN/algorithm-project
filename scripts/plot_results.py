import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    plt.style.use('ggplot')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--outdir', default='figs')
    ap.add_argument('--groupby', default='gen')
    args = ap.parse_args()

    df = pd.read_csv(args.csv, index_col=False)
    # Reset index if algo column was used as index
    if df.index.name == 'algo' or (df.columns[0] != 'algo' and 'algo' not in df.columns):
        df = df.reset_index()
        if df.columns[0] == 'index':
            df = df.rename(columns={'index': 'algo'})
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Separate Dijkstra baseline and delta-stepping results
    dijkstra_df = df[df['algo'] == 'dijkstra']
    delta_df = df[df['algo'] == 'delta']

    # Strong scaling per group (fixed input, vary threads)
    for key, g in delta_df.groupby(args.groupby):
        if 'threads' not in g.columns:
            continue
        
        # Get Dijkstra baseline time for this group
        dijk_base = dijkstra_df[dijkstra_df[args.groupby] == key]
        if not dijk_base.empty:
            dijkstra_time = dijk_base['time_ms'].iloc[0]
        else:
            # Fallback: use single-thread delta as baseline
            base = g[g['threads'] == 1]
            if base.empty:
                continue
            dijkstra_time = base['time_ms'].iloc[0]

        # Compute speedup relative to Dijkstra baseline
        g = g.copy()
        g['speedup_vs_dijkstra'] = dijkstra_time / g['time_ms']
        
        # Also compute speedup relative to single-thread delta (parallel efficiency)
        base_delta = g[g['threads'] == 1]
        if not base_delta.empty:
            for d in g['delta'].unique():
                mask = g['delta'] == d
                base_time = base_delta[base_delta['delta'] == d]['time_ms'].values
                if len(base_time) > 0:
                    g.loc[mask, 'speedup_vs_delta1'] = base_time[0] / g.loc[mask, 'time_ms']
        
        g['efficiency'] = g.get('speedup_vs_delta1', g['speedup_vs_dijkstra']) / g['threads']

        # Plot speedup vs Dijkstra for each delta
        for d, dfg in g.groupby('delta'):
            dfg = dfg.sort_values('threads')
            
            # Speedup vs Dijkstra
            plt.figure(figsize=(8, 5))
            plt.plot(dfg['threads'], dfg['speedup_vs_dijkstra'], marker='o', linewidth=2, markersize=8)
            plt.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Dijkstra baseline')
            plt.xlabel('Threads', fontsize=12)
            plt.ylabel('Speedup (vs Dijkstra)', fontsize=12)
            plt.title(f'Strong Scaling: {key.upper()}, Δ={d}', fontsize=14)
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(Path(args.outdir) / f'speedup_{key}_delta{d}.png', bbox_inches='tight', dpi=150)
            plt.close()

            # Parallel efficiency
            plt.figure(figsize=(8, 5))
            plt.plot(dfg['threads'], dfg['efficiency'], marker='s', linewidth=2, markersize=8, color='green')
            plt.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Ideal efficiency')
            plt.xlabel('Threads', fontsize=12)
            plt.ylabel('Efficiency (Sp/P)', fontsize=12)
            plt.ylim(0, 1.2)
            plt.title(f'Parallel Efficiency: {key.upper()}, Δ={d}', fontsize=14)
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(Path(args.outdir) / f'efficiency_{key}_delta{d}.png', bbox_inches='tight', dpi=150)
            plt.close()

        # Delta sensitivity plot (time vs delta for different thread counts)
        plt.figure(figsize=(8, 5))
        for t in sorted(g['threads'].unique()):
            subset = g[g['threads'] == t].sort_values('delta')
            plt.plot(subset['delta'], subset['time_ms'], marker='o', label=f'{t} threads')
        plt.xlabel('Delta (Δ)', fontsize=12)
        plt.ylabel('Time (ms)', fontsize=12)
        plt.title(f'Delta Sensitivity: {key.upper()}', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(Path(args.outdir) / f'delta_sensitivity_{key}.png', bbox_inches='tight', dpi=150)
        plt.close()

    # Summary comparison across graph types
    if len(delta_df[args.groupby].unique()) > 1:
        plt.figure(figsize=(10, 6))
        for key in delta_df[args.groupby].unique():
            subset = delta_df[(delta_df[args.groupby] == key) & (delta_df['delta'] == delta_df['delta'].median())]
            if not subset.empty:
                subset = subset.sort_values('threads')
                plt.plot(subset['threads'], subset['time_ms'], marker='o', label=key.upper())
        plt.xlabel('Threads', fontsize=12)
        plt.ylabel('Time (ms)', fontsize=12)
        plt.title('Runtime Comparison Across Graph Types', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(Path(args.outdir) / 'comparison_all.png', bbox_inches='tight', dpi=150)
        plt.close()

    print(f"Plots saved to {args.outdir}/")

if __name__ == '__main__':
    main()
