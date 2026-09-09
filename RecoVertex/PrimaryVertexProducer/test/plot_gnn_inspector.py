#!/usr/bin/env python3
"""
plot_gnn_inspector.py - Draw all histograms from GNN Track Inspector output
"""

import sys
import uproot
import matplotlib.pyplot as plt
import numpy as np

def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'v17s7'
    # Open file
    f = uproot.open(f'gnn_inspector_{tag}.root')
    
    # Get all histogram keys (skip directory)
    hist_keys = [k for k in f.keys() if '/' in k]
    
    print(f"Found {len(hist_keys)} histograms:")
    for k in hist_keys:
        print(f"  {k}")
    
    # Create figure with subplots
    n_hists = len(hist_keys)
    n_cols = 3
    n_rows = (n_hists + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    for i, key in enumerate(hist_keys):
        h = f[key]
        
        # Get histogram name
        hist_name = key.split('/')[-1].replace(';1', '')
        
        # Get values and edges
        values = h.values()
        edges = h.axis().edges()
        centers = (edges[:-1] + edges[1:]) / 2
        
        # Plot
        ax = axes[i]
        ax.bar(centers, values, width=np.diff(edges), alpha=0.7, edgecolor='black')
        ax.set_title(hist_name, fontsize=14, fontweight='bold')
        ax.set_xlabel('Value')
        ax.set_ylabel('Entries')
        ax.grid(True, alpha=0.3)
        
        # Add stats
        total = np.sum(values)
        mean = np.average(centers, weights=values) if total > 0 else 0
        ax.text(0.95, 0.95, f'Entries: {int(total)}\nMean: {mean:.3f}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Hide empty subplots
    for i in range(len(hist_keys), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'gnn_inspector_histograms_{tag}.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved to gnn_inspector_histograms_{tag}.png")
    plt.close()

if __name__ == "__main__":
    main()
