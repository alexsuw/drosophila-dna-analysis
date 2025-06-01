#!/usr/bin/env python3
"""
Z-Hunt Results Analysis Script
Analyzes Z-DNA predictions from Z-Hunt and filters by Z-score threshold
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import argparse
import os

def parse_zhunt_output(zhunt_file):
    """
    Parse Z-Hunt output file and extract relevant information
    """
    results = []
    
    with open(zhunt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Z-Hunt output format: chr start end length score sequence
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        result = {
                            'chromosome': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'length': int(parts[3]),
                            'z_score': float(parts[4]),
                            'sequence': parts[5] if len(parts) > 5 else ''
                        }
                        results.append(result)
                    except (ValueError, IndexError):
                        continue
    
    return pd.DataFrame(results)

def filter_by_zscore(df, threshold_min=300, threshold_max=400):
    """
    Filter Z-DNA results by Z-score threshold
    """
    filtered_df = df[(df['z_score'] >= threshold_min) & (df['z_score'] <= threshold_max)]
    print(f"Filtered {len(filtered_df)} regions with Z-score between {threshold_min} and {threshold_max}")
    print(f"Original: {len(df)} regions, Filtered: {len(filtered_df)} regions")
    return filtered_df

def analyze_zdna_distribution(df):
    """
    Analyze Z-DNA distribution across chromosomes
    """
    # Basic statistics
    print("\n=== Z-DNA Distribution Analysis ===")
    print(f"Total Z-DNA regions: {len(df)}")
    print(f"Average Z-score: {df['z_score'].mean():.2f}")
    print(f"Median Z-score: {df['z_score'].median():.2f}")
    print(f"Z-score range: {df['z_score'].min():.2f} - {df['z_score'].max():.2f}")
    
    # Distribution by chromosome
    chr_counts = df['chromosome'].value_counts()
    print(f"\nZ-DNA regions per chromosome:")
    print(chr_counts.head(10))
    
    return chr_counts

def create_bed_file(df, output_file):
    """
    Create BED format file for further analysis with bedtools
    """
    bed_df = df[['chromosome', 'start', 'end', 'z_score']].copy()
    bed_df['name'] = 'Z-DNA_' + bed_df.index.astype(str)
    bed_df['strand'] = '.'
    
    # BED format: chr start end name score strand
    bed_df = bed_df[['chromosome', 'start', 'end', 'name', 'z_score', 'strand']]
    bed_df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"BED file saved to: {output_file}")

def plot_zscore_distribution(df, output_dir):
    """
    Create plots for Z-score distribution analysis
    """
    plt.figure(figsize=(15, 10))
    
    # Z-score histogram
    plt.subplot(2, 3, 1)
    plt.hist(df['z_score'], bins=50, alpha=0.7, color='skyblue')
    plt.xlabel('Z-score')
    plt.ylabel('Frequency')
    plt.title('Z-score Distribution')
    plt.grid(True, alpha=0.3)
    
    # Length distribution
    plt.subplot(2, 3, 2)
    plt.hist(df['length'], bins=30, alpha=0.7, color='lightgreen')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Z-DNA Region Length Distribution')
    plt.grid(True, alpha=0.3)
    
    # Z-score vs Length scatter plot
    plt.subplot(2, 3, 3)
    plt.scatter(df['length'], df['z_score'], alpha=0.6, s=10)
    plt.xlabel('Length (bp)')
    plt.ylabel('Z-score')
    plt.title('Z-score vs Length')
    plt.grid(True, alpha=0.3)
    
    # Chromosome distribution
    plt.subplot(2, 3, 4)
    chr_counts = df['chromosome'].value_counts().head(10)
    chr_counts.plot(kind='bar')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Z-DNA regions')
    plt.title('Z-DNA Distribution by Chromosome')
    plt.xticks(rotation=45)
    
    # Z-score box plot by chromosome
    plt.subplot(2, 3, 5)
    top_chrs = df['chromosome'].value_counts().head(8).index
    df_top_chrs = df[df['chromosome'].isin(top_chrs)]
    sns.boxplot(data=df_top_chrs, x='chromosome', y='z_score')
    plt.xlabel('Chromosome')
    plt.ylabel('Z-score')
    plt.title('Z-score Distribution by Chromosome')
    plt.xticks(rotation=45)
    
    # Cumulative Z-score distribution
    plt.subplot(2, 3, 6)
    sorted_scores = np.sort(df['z_score'])
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    plt.plot(sorted_scores, cumulative)
    plt.xlabel('Z-score')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Z-score Distribution')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/zdna_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Analysis plots saved to: {output_dir}/zdna_analysis.png")

def main():
    parser = argparse.ArgumentParser(description='Analyze Z-Hunt results')
    parser.add_argument('--input', required=True, help='Z-Hunt output file')
    parser.add_argument('--output-dir', default='results', help='Output directory')
    parser.add_argument('--min-zscore', type=float, default=300, help='Minimum Z-score threshold')
    parser.add_argument('--max-zscore', type=float, default=400, help='Maximum Z-score threshold')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=== Z-Hunt Results Analysis ===")
    print(f"Input file: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"Z-score filter: {args.min_zscore} - {args.max_zscore}")
    
    # Parse Z-Hunt results
    print("\nParsing Z-Hunt output...")
    df = parse_zhunt_output(args.input)
    
    if df.empty:
        print("No Z-DNA regions found in the input file!")
        return
    
    # Analyze distribution before filtering
    print("\n=== Before Filtering ===")
    analyze_zdna_distribution(df)
    
    # Filter by Z-score
    print(f"\n=== Filtering by Z-score ({args.min_zscore}-{args.max_zscore}) ===")
    filtered_df = filter_by_zscore(df, args.min_zscore, args.max_zscore)
    
    if filtered_df.empty:
        print("No Z-DNA regions pass the Z-score filter!")
        return
    
    # Analyze filtered results
    print("\n=== After Filtering ===")
    analyze_zdna_distribution(filtered_df)
    
    # Save results
    filtered_df.to_csv(f"{args.output_dir}/zdna_filtered.csv", index=False)
    create_bed_file(filtered_df, f"{args.output_dir}/zdna_filtered.bed")
    
    # Create visualizations
    print("\nCreating visualizations...")
    plot_zscore_distribution(filtered_df, args.output_dir)
    
    print(f"\nAnalysis complete! Check {args.output_dir}/ for results.")

if __name__ == "__main__":
    main() 