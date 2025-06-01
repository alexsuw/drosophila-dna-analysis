#!/usr/bin/env python3
"""
G-Quadruplex Search Script
Searches for G-quadruplex forming sequences in genome using pattern matching
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
from collections import defaultdict

def find_gquadruplex_patterns(sequence, min_run_length=3, max_loop_length=7):
    """
    Find G-quadruplex forming sequences using regex patterns
    G-quadruplex pattern: G{n}N{1-7}G{n}N{1-7}G{n}N{1-7}G{n}
    where n >= 3 (minimum run length) and N is any nucleotide
    """
    patterns = []
    
    # Generate patterns for different G-run lengths (3-6)
    for g_run in range(min_run_length, 7):
        # Pattern: G{n}N{1-max_loop}G{n}N{1-max_loop}G{n}N{1-max_loop}G{n}
        pattern = f"G{{{g_run},}}[ATCGN]{{1,{max_loop_length}}}G{{{g_run},}}[ATCGN]{{1,{max_loop_length}}}G{{{g_run},}}[ATCGN]{{1,{max_loop_length}}}G{{{g_run},}}"
        patterns.append((pattern, g_run))
    
    return patterns

def search_quadruplexes_in_sequence(seq_record, patterns, chromosome_name):
    """
    Search for G-quadruplex patterns in a single sequence
    """
    results = []
    sequence = str(seq_record.seq).upper()
    
    for pattern_str, g_run_length in patterns:
        pattern = re.compile(pattern_str)
        
        for match in pattern.finditer(sequence):
            start = match.start()
            end = match.end()
            matched_seq = match.group()
            
            # Calculate some basic properties
            g_content = matched_seq.count('G') / len(matched_seq)
            gc_content = (matched_seq.count('G') + matched_seq.count('C')) / len(matched_seq)
            
            result = {
                'chromosome': chromosome_name,
                'start': start,
                'end': end,
                'length': end - start,
                'sequence': matched_seq,
                'g_run_length': g_run_length,
                'g_content': g_content,
                'gc_content': gc_content,
                'score': calculate_gquad_score(matched_seq)
            }
            results.append(result)
    
    return results

def calculate_gquad_score(sequence):
    """
    Calculate a simple scoring for G-quadruplex potential
    Based on G-content, G-runs, and sequence properties
    """
    score = 0
    
    # G-content contribution (higher G-content = higher score)
    g_content = sequence.count('G') / len(sequence)
    score += g_content * 100
    
    # Count G-runs and their lengths
    g_runs = re.findall(r'G+', sequence)
    if len(g_runs) >= 4:  # At least 4 G-runs for quadruplex
        score += len(g_runs) * 10
        # Bonus for longer G-runs
        avg_run_length = sum(len(run) for run in g_runs) / len(g_runs)
        score += avg_run_length * 5
    
    # Penalty for very long loops (less stable quadruplexes)
    non_g_stretches = re.findall(r'[ATC]+', sequence)
    if non_g_stretches:
        avg_loop_length = sum(len(loop) for loop in non_g_stretches) / len(non_g_stretches)
        if avg_loop_length > 5:
            score *= 0.8  # Reduce score for long loops
    
    return round(score, 2)

def search_genome_quadruplexes(fasta_file, min_score=50):
    """
    Search for G-quadruplexes in the entire genome
    """
    print(f"Searching for G-quadruplexes in {fasta_file}")
    
    patterns = find_gquadruplex_patterns("", min_run_length=3, max_loop_length=7)
    all_results = []
    
    # Parse FASTA file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = seq_record.id
        print(f"Processing chromosome: {chromosome}")
        
        # Search for quadruplexes
        results = search_quadruplexes_in_sequence(seq_record, patterns, chromosome)
        all_results.extend(results)
        
        print(f"  Found {len(results)} potential G-quadruplexes")
    
    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    
    if df.empty:
        print("No G-quadruplexes found!")
        return df
    
    # Filter by score
    filtered_df = df[df['score'] >= min_score]
    print(f"\nFiltered {len(filtered_df)} G-quadruplexes with score >= {min_score}")
    print(f"Original: {len(df)}, Filtered: {len(filtered_df)}")
    
    return filtered_df

def analyze_quadruplex_distribution(df):
    """
    Analyze G-quadruplex distribution and properties
    """
    print("\n=== G-Quadruplex Distribution Analysis ===")
    print(f"Total G-quadruplex regions: {len(df)}")
    print(f"Average score: {df['score'].mean():.2f}")
    print(f"Median score: {df['score'].median():.2f}")
    print(f"Score range: {df['score'].min():.2f} - {df['score'].max():.2f}")
    
    print(f"\nAverage length: {df['length'].mean():.1f} bp")
    print(f"Length range: {df['length'].min()} - {df['length'].max()} bp")
    
    print(f"\nAverage G-content: {df['g_content'].mean():.3f}")
    print(f"Average GC-content: {df['gc_content'].mean():.3f}")
    
    # Distribution by chromosome
    chr_counts = df['chromosome'].value_counts()
    print(f"\nG-quadruplex regions per chromosome:")
    print(chr_counts.head(10))
    
    return chr_counts

def create_bed_file(df, output_file):
    """
    Create BED format file for G-quadruplexes
    """
    bed_df = df[['chromosome', 'start', 'end', 'score']].copy()
    bed_df['name'] = 'G4_' + bed_df.index.astype(str)
    bed_df['strand'] = '.'
    
    # BED format: chr start end name score strand
    bed_df = bed_df[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    bed_df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"BED file saved to: {output_file}")

def plot_quadruplex_analysis(df, output_dir):
    """
    Create plots for G-quadruplex analysis
    """
    plt.figure(figsize=(15, 10))
    
    # Score distribution
    plt.subplot(2, 3, 1)
    plt.hist(df['score'], bins=50, alpha=0.7, color='orange')
    plt.xlabel('G4 Score')
    plt.ylabel('Frequency')
    plt.title('G-Quadruplex Score Distribution')
    plt.grid(True, alpha=0.3)
    
    # Length distribution
    plt.subplot(2, 3, 2)
    plt.hist(df['length'], bins=30, alpha=0.7, color='lightcoral')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    plt.title('G-Quadruplex Length Distribution')
    plt.grid(True, alpha=0.3)
    
    # G-content vs Score
    plt.subplot(2, 3, 3)
    plt.scatter(df['g_content'], df['score'], alpha=0.6, s=10)
    plt.xlabel('G-content')
    plt.ylabel('G4 Score')
    plt.title('G-content vs Score')
    plt.grid(True, alpha=0.3)
    
    # Chromosome distribution
    plt.subplot(2, 3, 4)
    chr_counts = df['chromosome'].value_counts().head(10)
    chr_counts.plot(kind='bar')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of G-quadruplexes')
    plt.title('G-Quadruplex Distribution by Chromosome')
    plt.xticks(rotation=45)
    
    # G-run length distribution
    plt.subplot(2, 3, 5)
    plt.hist(df['g_run_length'], bins=range(3, 8), alpha=0.7, color='lightblue')
    plt.xlabel('G-run Length')
    plt.ylabel('Frequency')
    plt.title('G-run Length Distribution')
    plt.grid(True, alpha=0.3)
    
    # GC content distribution
    plt.subplot(2, 3, 6)
    plt.hist(df['gc_content'], bins=30, alpha=0.7, color='lightgreen')
    plt.xlabel('GC Content')
    plt.ylabel('Frequency')
    plt.title('GC Content Distribution')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/quadruplex_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Analysis plots saved to: {output_dir}/quadruplex_analysis.png")

def main():
    parser = argparse.ArgumentParser(description='Search for G-quadruplexes in genome')
    parser.add_argument('--input', required=True, help='Genome FASTA file')
    parser.add_argument('--output-dir', default='results', help='Output directory')
    parser.add_argument('--min-score', type=float, default=50, help='Minimum G4 score threshold')
    parser.add_argument('--min-run-length', type=int, default=3, help='Minimum G-run length')
    parser.add_argument('--max-loop-length', type=int, default=7, help='Maximum loop length')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=== G-Quadruplex Search ===")
    print(f"Input file: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum score: {args.min_score}")
    print(f"G-run length: >= {args.min_run_length}")
    print(f"Max loop length: {args.max_loop_length}")
    
    # Search for G-quadruplexes
    df = search_genome_quadruplexes(args.input, args.min_score)
    
    if df.empty:
        print("No G-quadruplexes found!")
        return
    
    # Analyze results
    analyze_quadruplex_distribution(df)
    
    # Save results
    df.to_csv(f"{args.output_dir}/quadruplex_results.csv", index=False)
    create_bed_file(df, f"{args.output_dir}/quadruplex_results.bed")
    
    # Create visualizations
    print("\nCreating visualizations...")
    plot_quadruplex_analysis(df, args.output_dir)
    
    print(f"\nG-quadruplex search complete! Check {args.output_dir}/ for results.")

if __name__ == "__main__":
    main() 