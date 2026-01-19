#!/usr/bin/env python3
"""
TC2038 - Benchmark Visualization Script
Generates plots for algorithm runtime analysis
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def plot_string_benchmarks():
    """Plot string algorithm benchmarks."""
    if not os.path.exists('string_benchmark.csv'):
        print("Warning: string_benchmark.csv not found. Run benchmark first.")
        return
    
    df = pd.read_csv('string_benchmark.csv')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Linear scale plot
    ax1.plot(df['size'], df['kmp_ms'], 'o-', label='KMP', linewidth=2, markersize=8)
    ax1.plot(df['size'], df['z_algo_ms'], 's-', label='Z-Algorithm', linewidth=2, markersize=8)
    ax1.plot(df['size'], df['rabin_karp_ms'], '^-', label='Rabin-Karp', linewidth=2, markersize=8)
    ax1.set_xlabel('Text Size (characters)', fontsize=12)
    ax1.set_ylabel('Time (ms)', fontsize=12)
    ax1.set_title('String Matching Algorithms - Runtime Comparison', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Log-log scale plot (to verify O(n) complexity)
    ax2.loglog(df['size'], df['kmp_ms'], 'o-', label='KMP', linewidth=2, markersize=8)
    ax2.loglog(df['size'], df['z_algo_ms'], 's-', label='Z-Algorithm', linewidth=2, markersize=8)
    ax2.loglog(df['size'], df['rabin_karp_ms'], '^-', label='Rabin-Karp', linewidth=2, markersize=8)
    # Reference line for O(n)
    sizes = df['size'].values
    ref_line = sizes * (df['kmp_ms'].iloc[0] / sizes[0])
    ax2.loglog(sizes, ref_line, '--', color='gray', alpha=0.7, label='O(n) reference')
    ax2.set_xlabel('Text Size (characters)', fontsize=12)
    ax2.set_ylabel('Time (ms)', fontsize=12)
    ax2.set_title('Log-Log Scale (slope ≈ 1 confirms O(n))', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('string_benchmark_plot.png', dpi=150, bbox_inches='tight')
    print("Saved: string_benchmark_plot.png")
    plt.close()

def plot_dna_benchmarks():
    """Plot DNA algorithm benchmarks."""
    if not os.path.exists('dna_benchmark.csv'):
        print("Warning: dna_benchmark.csv not found. Run benchmark first.")
        return
    
    df = pd.read_csv('dna_benchmark.csv')
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Manacher
    ax1 = axes[0]
    ax1.plot(df['size'], df['manacher_ms'], 'o-', color='tab:blue', linewidth=2, markersize=8)
    ax1.set_xlabel('Sequence Size', fontsize=12)
    ax1.set_ylabel('Time (ms)', fontsize=12)
    ax1.set_title('Manacher\'s Algorithm - O(n)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Suffix Array Construction
    ax2 = axes[1]
    ax2.plot(df['size'], df['suffix_array_ms'], 's-', color='tab:orange', linewidth=2, markersize=8)
    ax2.set_xlabel('Sequence Size', fontsize=12)
    ax2.set_ylabel('Time (ms)', fontsize=12)
    ax2.set_title('Suffix Array Construction - O(n log n)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # LCS DP
    ax3 = axes[2]
    ax3.plot(df['size'], df['lcs_dp_ms'], 'd-', color='tab:red', linewidth=2, markersize=8)
    ax3.set_xlabel('Sequence Size', fontsize=12)
    ax3.set_ylabel('Time (ms)', fontsize=12)
    ax3.set_title('LCS via Dynamic Programming - O(nm)', fontsize=14)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('dna_benchmark_plot.png', dpi=150, bbox_inches='tight')
    print("Saved: dna_benchmark_plot.png")
    plt.close()

def plot_graph_benchmarks():
    """Plot graph algorithm benchmarks."""
    if not os.path.exists('graph_benchmark.csv'):
        print("Warning: graph_benchmark.csv not found. Run benchmark first.")
        return
    
    df = pd.read_csv('graph_benchmark.csv')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # MST comparison - Linear scale
    ax1.plot(df['vertices'], df['kruskal_ms'], '^-', label='Kruskal', linewidth=2, markersize=8)
    ax1.plot(df['vertices'], df['prim_ms'], 'd-', label='Prim', linewidth=2, markersize=8)
    ax1.set_xlabel('Number of Vertices', fontsize=12)
    ax1.set_ylabel('Time (ms)', fontsize=12)
    ax1.set_title('MST Algorithms - Runtime Comparison', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # MST comparison - Log-log scale
    ax2.loglog(df['vertices'], df['kruskal_ms'], '^-', label='Kruskal', linewidth=2, markersize=8)
    ax2.loglog(df['vertices'], df['prim_ms'], 'd-', label='Prim', linewidth=2, markersize=8)
    # Reference line for O(E log E) where E = 3*V
    vertices = df['vertices'].values
    ref_kruskal = vertices * np.log(3*vertices) * (df['kruskal_ms'].iloc[0] / (vertices[0] * np.log(3*vertices[0])))
    ax2.loglog(vertices, ref_kruskal, '--', color='gray', alpha=0.7, label='O(E log E) reference')
    ax2.set_xlabel('Number of Vertices', fontsize=12)
    ax2.set_ylabel('Time (ms)', fontsize=12)
    ax2.set_title('Log-Log Scale (verifies complexity)', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('graph_benchmark_plot.png', dpi=150, bbox_inches='tight')
    print("Saved: graph_benchmark_plot.png")
    plt.close()

def create_summary_table():
    """Create a summary table of algorithm complexities."""
    data = {
        'Algorithm': [
            'KMP', 'Z-Algorithm', 'Rabin-Karp',
            'Manacher', 'Suffix Array + LCP', 'LCS (SA)', 'LCS (DP)',
            'Edmonds-Karp', 'Dinic', 'Kruskal', 'Prim',
            'Second MST', 'Biconnected', 'Expression Parser'
        ],
        'Time Complexity': [
            'O(n + m)', 'O(n + m)', 'O(n + m) avg',
            'O(n)', 'O(n log n)', 'O(n log n)', 'O(nm)',
            'O(VE²)', 'O(V²E)', 'O(E log E)', 'O(E log V)',
            'O(E log E + V²)', 'O(V + E)', 'O(n)'
        ],
        'Space Complexity': [
            'O(m)', 'O(n)', 'O(m)',
            'O(n)', 'O(n)', 'O(n)', 'O(nm)',
            'O(V²)', 'O(V²)', 'O(V)', 'O(V)',
            'O(V + E)', 'O(V + E)', 'O(n)'
        ],
        'Problem': [
            'String matching', 'String matching', 'String matching',
            'Palindromes', 'Suffix/LCP', 'LCS', 'LCS',
            'Max Flow', 'Max Flow', 'MST', 'MST',
            '2nd MST', 'Bridges/BCC', 'Expression eval'
        ]
    }
    
    df = pd.DataFrame(data)
    df.to_csv('complexity_summary.csv', index=False)
    print("Saved: complexity_summary.csv")
    
    # Print as table
    print("\n" + "="*80)
    print("ALGORITHM COMPLEXITY SUMMARY")
    print("="*80)
    print(df.to_string(index=False))
    print("="*80)

def main():
    print("TC2038 - Benchmark Visualization")
    print("-" * 40)
    
    # Check for required libraries
    try:
        import matplotlib
        import pandas
    except ImportError as e:
        print(f"Error: {e}")
        print("Please install required packages: pip install matplotlib pandas")
        return
    
    # Generate all plots
    plot_string_benchmarks()
    plot_dna_benchmarks()
    plot_graph_benchmarks()
    create_summary_table()
    
    print("\nDone! Check the generated PNG files for plots.")
    print("Files generated:")
    for f in ['string_benchmark_plot.png', 'dna_benchmark_plot.png', 
              'graph_benchmark_plot.png', 'complexity_summary.csv']:
        if os.path.exists(f):
            print(f"  ✓ {f}")
        else:
            print(f"  ✗ {f} (benchmark CSV missing)")

if __name__ == "__main__":
    main()
