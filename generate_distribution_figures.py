#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate 6 high-quality distribution figures from ATLAS database.

Figures:
(b) Motif Type Distribution
(c) Internal Loop Size Distribution
(d) Hairpin Loop Size Distribution
(e) Bulge Size Distribution
(f) 3-way Junction Size Distribution
(g) Pseudoknot Type Distribution
"""

import sqlite3
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# Database path
base_dir = os.path.dirname(__file__)
db_path = os.path.join(base_dir, 'ATLAS.db')

# Global style settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5

def apply_global_style(ax, ylabel="Number of Motifs"):
    """Apply consistent styling to all figures."""
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add faint dashed horizontal grid
    ax.grid(axis='y', alpha=0.5, linestyle='--', linewidth=0.8)
    ax.set_axisbelow(True)

    # Set y-axis label with larger font
    ax.set_ylabel(ylabel, fontsize=24, fontweight='bold')

    # Increase tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=21)

    # Thicken remaining spines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

def annotate_bars(ax, bars, fontsize=9):
    """Add value labels on top of bars."""
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height):,}',
                   ha='center', va='bottom', fontsize=fontsize, fontweight='bold')

def plot_figure_b(df):
    """(b) Motif Type Distribution - Categorical Bar Chart"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Count motif types
    motif_counts = df['motif_type'].value_counts().sort_values(ascending=False)

    # Custom gradient colors (Dark Teal -> Teal -> Orange/Yellow)
    colors = ['#264653', '#2A9D8F', '#2A9D8F', '#E76F51', '#F4A261',
              '#E9C46A', '#E9C46A', '#E9C46A', '#E9C46A', '#E9C46A', '#E9C46A']
    colors = colors[:len(motif_counts)]

    # Create bar chart
    bars = ax.bar(range(len(motif_counts)), motif_counts.values, color=colors,
                  edgecolor='white', linewidth=1.5)

    # Set x-axis labels with larger font and rotation
    ax.set_xticks(range(len(motif_counts)))
    ax.set_xticklabels(motif_counts.index, rotation=45, ha='right', fontsize=21)
    ax.set_xlabel('Motif Type', fontsize=24, fontweight='bold')

    # Set y-axis limit - NO "k" suffix, use full numbers
    ax.set_ylim(0, 300000)
    ax.set_yticks(np.arange(0, 300001, 50000))

    # Apply global style
    apply_global_style(ax)

    # Annotate bars
    annotate_bars(ax, bars, fontsize=8)

    # No title
    plt.tight_layout()
    plt.savefig('figure_b_motif_type.png', dpi=300, bbox_inches='tight')
    print("Figure (b) saved: figure_b_motif_type.png")
    plt.close()

def plot_figure_c(df):
    """(c) Internal Loop Size Distribution - Histogram"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter internal loops and calculate LOOP size (not total)
    internal_df = df[df['motif_type'] == 'internal'].copy()
    internal_df['loop_size'] = internal_df['nt_number'].apply(
        lambda x: len(x.split(',')) - 4 if pd.notna(x) else 0  # Subtract 4 paired nucleotides
    )

    # Bin sizes 0-20
    sizes = internal_df['loop_size'].values
    sizes = sizes[(sizes >= 0) & (sizes <= 20)]

    # Create histogram
    bins = np.arange(0, 21, 1)
    counts, edges = np.histogram(sizes, bins=bins)

    bars = ax.bar(edges[:-1], counts, width=0.8, color='#2A9D8F',
                  edgecolor='white', linewidth=1.2)

    # Set axes - Y max changed to 100,000
    ax.set_xlim(-0.5, 20.5)
    ax.set_ylim(0, 100000)
    ax.set_xlabel('Internal Loop Size', fontsize=24, fontweight='bold')
    ax.set_xticks(range(0, 21, 2))
    ax.set_yticks(np.arange(0, 100001, 20000))

    # Apply global style
    apply_global_style(ax)

    # No annotations for Internal (font too small to read)

    # No title
    plt.tight_layout()
    plt.savefig('figure_c_internal_size.png', dpi=300, bbox_inches='tight')
    print("Figure (c) saved: figure_c_internal_size.png")
    plt.close()

def plot_figure_d(df):
    """(d) Hairpin Loop Size Distribution - Histogram"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter hairpin loops and calculate LOOP size (not total)
    hairpin_df = df[df['motif_type'] == 'hairpin'].copy()
    hairpin_df['loop_size'] = hairpin_df['nt_number'].apply(
        lambda x: len(x.split(',')) - 2 if pd.notna(x) else 0  # Subtract 2 paired nucleotides
    )

    # Bin sizes 0-20
    sizes = hairpin_df['loop_size'].values
    sizes = sizes[(sizes >= 0) & (sizes <= 20)]

    # Create histogram
    bins = np.arange(0, 21, 1)
    counts, edges = np.histogram(sizes, bins=bins)

    bars = ax.bar(edges[:-1], counts, width=0.8, color='#F4A261',
                  edgecolor='white', linewidth=1.2)

    # Set axes
    ax.set_xlim(-0.5, 20.5)
    ax.set_ylim(0, 35000)
    ax.set_xlabel('Hairpin Loop Size', fontsize=24, fontweight='bold')
    ax.set_xticks(range(0, 21, 2))
    ax.set_yticks(np.arange(0, 35001, 5000))

    # Apply global style
    apply_global_style(ax)

    # No annotations for Hairpin (font too small to read)

    # No title
    plt.tight_layout()
    plt.savefig('figure_d_hairpin_size.png', dpi=300, bbox_inches='tight')
    print("Figure (d) saved: figure_d_hairpin_size.png")
    plt.close()

def plot_figure_e(df):
    """(e) Bulge Size Distribution - Histogram"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter bulge loops and calculate LOOP size (not total)
    bulge_df = df[df['motif_type'] == 'bulge'].copy()
    bulge_df['loop_size'] = bulge_df['nt_number'].apply(
        lambda x: len(x.split(',')) - 4 if pd.notna(x) else 0  # Subtract 4 paired nucleotides
    )

    # Bin sizes 0-10 (extended range to check data)
    sizes = bulge_df['loop_size'].values
    sizes = sizes[(sizes >= 0) & (sizes <= 10)]

    # Create histogram
    bins = np.arange(0, 11, 1)
    counts, edges = np.histogram(sizes, bins=bins)

    bars = ax.bar(edges[:-1], counts, width=0.6, color='#E9C46A',
                  edgecolor='white', linewidth=1.5)

    # Set axes
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(0, max(counts) * 1.15 if max(counts) > 0 else 1000)
    ax.set_xlabel('Bulge Size', fontsize=24, fontweight='bold')
    ax.set_xticks(range(0, 11))

    # Apply global style
    apply_global_style(ax)

    # No annotations for Bulge (font too small to read)

    # No title
    plt.tight_layout()
    plt.savefig('figure_e_bulge_size.png', dpi=300, bbox_inches='tight')
    print("Figure (e) saved: figure_e_bulge_size.png")
    plt.close()

def plot_figure_f(df):
    """(f) 3-way Junction Size Distribution - Histogram"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter 3-way junctions and calculate LOOP size (not total)
    junction3_df = df[df['motif_type'] == '3_junction'].copy()
    junction3_df['loop_size'] = junction3_df['nt_number'].apply(
        lambda x: len(x.split(',')) - 6 if pd.notna(x) else 0  # Subtract 6 paired nucleotides
    )

    # Bin sizes
    sizes = junction3_df['loop_size'].values
    sizes = sizes[(sizes >= 0) & (sizes <= 25)]

    # Create histogram
    bins = np.arange(0, 26, 1)
    counts, edges = np.histogram(sizes, bins=bins)

    bars = ax.bar(edges[:-1], counts, width=0.8, color='#E76F51',
                  edgecolor='white', linewidth=1.2)

    # Set axes
    ax.set_xlim(-0.5, 25.5)
    ax.set_ylim(0, max(counts) * 1.15 if max(counts) > 0 else 1000)
    ax.set_xlabel('3-way Junction Size', fontsize=24, fontweight='bold')
    ax.set_xticks(range(0, 26, 2))

    # Apply global style
    apply_global_style(ax)

    # No annotations for 3-way Junction (font too small to read)

    # No title
    plt.tight_layout()
    plt.savefig('figure_f_3way_size.png', dpi=300, bbox_inches='tight')
    print("Figure (f) saved: figure_f_3way_size.png")
    plt.close()

def plot_figure_g(df_pk):
    """(g) Pseudoknot Type Distribution - Categorical Bar Chart"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Count pseudoknot types
    pk_counts = df_pk['motif_type'].value_counts().sort_values(ascending=False)

    # Multi-color palette
    colors = ['#264653', '#2A9D8F', '#52B788', '#E9C46A', '#F4A261', '#E76F51']
    colors = colors[:len(pk_counts)]

    # Create bar chart
    bars = ax.bar(range(len(pk_counts)), pk_counts.values, color=colors,
                  edgecolor='white', linewidth=1.5)

    # Simplify labels: remove "pseudoknot_" prefix
    simplified_labels = [label.replace('pseudoknot_', '') for label in pk_counts.index]

    # Set x-axis labels - NO rotation, larger font (×1.5 = 21 → 31.5)
    ax.set_xticks(range(len(pk_counts)))
    ax.set_xticklabels(simplified_labels, rotation=0, ha='center', fontsize=32)
    ax.set_xlabel('Pseudoknot Type', fontsize=24, fontweight='bold')

    # Set y-axis limit - NO "k" suffix
    ax.set_ylim(0, 3500)
    ax.set_yticks(np.arange(0, 3501, 500))

    # Apply global style
    apply_global_style(ax)

    # Annotate bars with larger font (×2)
    annotate_bars(ax, bars, fontsize=18)

    # No title
    plt.tight_layout()
    plt.savefig('figure_g_pseudoknot_type.png', dpi=300, bbox_inches='tight')
    print("Figure (g) saved: figure_g_pseudoknot_type.png")
    plt.close()

def main():
    """Main function to generate all figures."""
    print("=" * 70)
    print("Generating Distribution Figures from ATLAS Database")
    print("=" * 70)

    # Connect to database
    print("\nConnecting to ATLAS.db...")
    conn = sqlite3.connect(db_path)

    # Load data table
    print("Loading data from 'data' table...")
    df = pd.read_sql_query("SELECT motif_type, nt_number FROM data", conn)
    print(f"  Loaded {len(df)} records from data table")

    # Load PK table
    print("Loading data from 'PK' table...")
    df_pk = pd.read_sql_query("SELECT motif_type FROM PK", conn)
    print(f"  Loaded {len(df_pk)} records from PK table")

    conn.close()

    print("\n" + "=" * 70)
    print("Generating figures...")
    print("=" * 70 + "\n")

    # Generate all figures
    plot_figure_b(df)
    plot_figure_c(df)
    plot_figure_d(df)
    plot_figure_e(df)
    plot_figure_f(df)
    plot_figure_g(df_pk)

    print("\n" + "=" * 70)
    print("All figures generated successfully!")
    print("=" * 70)
    print("\nOutput files:")
    print("  - figure_b_motif_type.png")
    print("  - figure_c_internal_size.png")
    print("  - figure_d_hairpin_size.png")
    print("  - figure_e_bulge_size.png")
    print("  - figure_f_3way_size.png")
    print("  - figure_g_pseudoknot_type.png")

if __name__ == '__main__':
    main()
