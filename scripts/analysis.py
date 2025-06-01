#!/usr/bin/env python3
"""
Genomic Location Analysis Script
Analyzes Z-DNA and G-quadruplex locations relative to genes and promoters
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import requests
import json
from collections import defaultdict
from pathlib import Path
import sys

def parse_gtf_file(gtf_file):
    """
    Parse GTF file to extract gene information
    """
    print(f"📚 Loading gene annotations from {gtf_file}...")
    
    genes = []
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type == 'transcript':
                chromosome = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # Parse attributes
                attributes = parts[8]
                gene_id = None
                gene_name = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                
                if gene_id:
                    genes.append({
                        'chromosome': chromosome,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_id': gene_id,
                        'gene_name': gene_name or gene_id,
                        'tss': start if strand == '+' else end
                    })
    
    return pd.DataFrame(genes)

def create_promoter_regions(genes_df, upstream=1000, downstream=1000):
    """
    Create promoter regions (±1000 bp from TSS)
    """
    print(f"🎯 Finding promoter overlaps (±{upstream}/{downstream}bp from TSS)...")
    
    promoters = []
    
    for _, gene in genes_df.iterrows():
        if gene['strand'] == '+':
            # TSS is at start position
            tss = gene['start']
            prom_start = max(1, tss - upstream)
            prom_end = tss + downstream
        else:
            # TSS is at end position
            tss = gene['end']
            prom_start = max(1, tss - downstream)
            prom_end = tss + upstream
        
        promoters.append({
            'chromosome': gene['chromosome'],
            'start': prom_start,
            'end': prom_end,
            'gene_id': gene['gene_id'],
            'gene_name': gene['gene_name'],
            'strand': gene['strand'],
            'tss': tss
        })
    
    return pd.DataFrame(promoters)

def find_overlaps(features_df, regions_df, feature_name="feature"):
    """
    Find overlaps between features and genomic regions
    """
    overlaps = []
    
    for _, feature in features_df.iterrows():
        for _, region in regions_df.iterrows():
            if (feature['chromosome'] == region['chromosome'] and
                not (feature['end'] < region['start'] or feature['start'] > region['end'])):
                
                overlap = {
                    'feature_chr': feature['chromosome'],
                    'feature_start': feature['start'],
                    'feature_end': feature['end'],
                    'region_chr': region['chromosome'],
                    'region_start': region['start'],
                    'region_end': region['end'],
                    'gene_id': region.get('gene_id', ''),
                    'gene_name': region.get('gene_name', ''),
                    'feature_type': feature_name
                }
                
                if 'z_score' in feature:
                    overlap['z_score'] = feature['z_score']
                if 'score' in feature:
                    overlap['score'] = feature['score']
                
                overlaps.append(overlap)
    
    return pd.DataFrame(overlaps)

def analyze_genomic_distribution(structures_df, genes_df, output_dir):
    """Analyze genomic distribution of structures"""
    print("📊 Analyzing genomic distribution...")
    
    # Chromosome distribution
    chrom_dist = structures_df['chromosome'].value_counts()
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Chromosome distribution
    top_chroms = chrom_dist.head(10)
    axes[0, 0].bar(range(len(top_chroms)), top_chroms.values)
    axes[0, 0].set_xticks(range(len(top_chroms)))
    axes[0, 0].set_xticklabels(top_chroms.index, rotation=45)
    axes[0, 0].set_title('Distribution across Chromosomes (Top 10)')
    axes[0, 0].set_ylabel('Count')
    
    # Length distribution
    if 'length' in structures_df.columns:
        axes[0, 1].hist(structures_df['length'], bins=50, alpha=0.7)
        axes[0, 1].set_title('Structure Length Distribution')
        axes[0, 1].set_xlabel('Length (bp)')
        axes[0, 1].set_ylabel('Count')
    
    # Score distribution
    if 'score' in structures_df.columns:
        axes[1, 0].hist(structures_df['score'], bins=50, alpha=0.7)
        axes[1, 0].set_title('Score Distribution')
        axes[1, 0].set_xlabel('Score')
        axes[1, 0].set_ylabel('Count')
    
    # GC content distribution
    if 'gc_content' in structures_df.columns:
        axes[1, 1].hist(structures_df['gc_content'], bins=50, alpha=0.7)
        axes[1, 1].set_title('GC Content Distribution')
        axes[1, 1].set_xlabel('GC Content')
        axes[1, 1].set_ylabel('Count')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'genomic_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return chrom_dist

def create_gene_lists(overlaps_dict):
    """
    Create gene lists for functional enrichment analysis
    """
    gene_lists = {}
    
    # Z-DNA genes
    zdna_genes = set()
    for _, overlap in overlaps_dict['zdna_gene_overlaps'].iterrows():
        if overlap['gene_id']:
            zdna_genes.add(overlap['gene_id'])
    
    # Z-DNA promoter genes
    zdna_promoter_genes = set()
    for _, overlap in overlaps_dict['zdna_promoter_overlaps'].iterrows():
        if overlap['gene_id']:
            zdna_promoter_genes.add(overlap['gene_id'])
    
    # G4 genes
    g4_genes = set()
    for _, overlap in overlaps_dict['g4_gene_overlaps'].iterrows():
        if overlap['gene_id']:
            g4_genes.add(overlap['gene_id'])
    
    # G4 promoter genes
    g4_promoter_genes = set()
    for _, overlap in overlaps_dict['g4_promoter_overlaps'].iterrows():
        if overlap['gene_id']:
            g4_promoter_genes.add(overlap['gene_id'])
    
    gene_lists = {
        'zdna_genes': list(zdna_genes),
        'zdna_promoter_genes': list(zdna_promoter_genes),
        'g4_genes': list(g4_genes),
        'g4_promoter_genes': list(g4_promoter_genes)
    }
    
    print(f"\n=== Gene Lists ===")
    print(f"Genes with Z-DNA: {len(gene_lists['zdna_genes'])}")
    print(f"Genes with Z-DNA in promoters: {len(gene_lists['zdna_promoter_genes'])}")
    print(f"Genes with G-quadruplexes: {len(gene_lists['g4_genes'])}")
    print(f"Genes with G-quadruplexes in promoters: {len(gene_lists['g4_promoter_genes'])}")
    
    return gene_lists

def save_gene_lists(gene_lists, output_dir):
    """
    Save gene lists for STRING DB analysis
    """
    for list_name, genes in gene_lists.items():
        output_file = f"{output_dir}/{list_name}.txt"
        with open(output_file, 'w') as f:
            for gene in genes:
                f.write(f"{gene}\n")
        print(f"Saved {len(genes)} genes to {output_file}")

def plot_genomic_analysis(overlaps_dict, gene_lists, output_dir):
    """
    Create visualization plots for genomic analysis
    """
    plt.figure(figsize=(15, 12))
    
    # Distribution pie charts
    plt.subplot(2, 3, 1)
    zdna_counts = [
        len(gene_lists['zdna_promoter_genes']),
        len(gene_lists['zdna_genes']) - len(gene_lists['zdna_promoter_genes']),
        # Approximate intergenic (this is simplified)
        max(0, len(overlaps_dict['zdna_gene_overlaps']) - len(gene_lists['zdna_genes']))
    ]
    labels = ['Promoters', 'Genes', 'Other']
    plt.pie(zdna_counts, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.title('Z-DNA Distribution')
    
    plt.subplot(2, 3, 2)
    g4_counts = [
        len(gene_lists['g4_promoter_genes']),
        len(gene_lists['g4_genes']) - len(gene_lists['g4_promoter_genes']),
        max(0, len(overlaps_dict['g4_gene_overlaps']) - len(gene_lists['g4_genes']))
    ]
    plt.pie(g4_counts, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.title('G-Quadruplex Distribution')
    
    # Gene count comparison
    plt.subplot(2, 3, 3)
    categories = ['Z-DNA\nGenes', 'Z-DNA\nPromoters', 'G4\nGenes', 'G4\nPromoters']
    counts = [
        len(gene_lists['zdna_genes']),
        len(gene_lists['zdna_promoter_genes']),
        len(gene_lists['g4_genes']),
        len(gene_lists['g4_promoter_genes'])
    ]
    bars = plt.bar(categories, counts, color=['skyblue', 'lightblue', 'orange', 'lightsalmon'])
    plt.ylabel('Number of Genes')
    plt.title('Gene Counts by Category')
    plt.xticks(rotation=45)
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                str(count), ha='center', va='bottom')
    
    # Overlap analysis
    plt.subplot(2, 3, 4)
    zdna_promoter_set = set(gene_lists['zdna_promoter_genes'])
    g4_promoter_set = set(gene_lists['g4_promoter_genes'])
    
    overlap_genes = zdna_promoter_set & g4_promoter_set
    zdna_only = zdna_promoter_set - g4_promoter_set
    g4_only = g4_promoter_set - zdna_promoter_set
    
    venn_data = [len(zdna_only), len(overlap_genes), len(g4_only)]
    venn_labels = ['Z-DNA only', 'Both', 'G4 only']
    plt.bar(venn_labels, venn_data, color=['skyblue', 'purple', 'orange'])
    plt.ylabel('Number of Genes')
    plt.title('Promoter Gene Overlap')
    plt.xticks(rotation=45)
    
    # Chromosome distribution
    plt.subplot(2, 3, 5)
    if not overlaps_dict['zdna_promoter_overlaps'].empty:
        zdna_chr_counts = overlaps_dict['zdna_promoter_overlaps']['feature_chr'].value_counts().head(8)
        zdna_chr_counts.plot(kind='bar', alpha=0.7, color='skyblue', label='Z-DNA')
    
    if not overlaps_dict['g4_promoter_overlaps'].empty:
        g4_chr_counts = overlaps_dict['g4_promoter_overlaps']['feature_chr'].value_counts().head(8)
        # Plot on same axes
        ax = plt.gca()
        g4_chr_counts.plot(kind='bar', alpha=0.7, color='orange', label='G4', ax=ax)
    
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Features')
    plt.title('Promoter Features by Chromosome')
    plt.legend()
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/genomic_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Genomic analysis plots saved to: {output_dir}/genomic_analysis.png")

def main():
    parser = argparse.ArgumentParser(description='Genomic analysis of Z-DNA and G-quadruplex structures')
    parser.add_argument('--zdna-file', help='Z-DNA results CSV file', default=None)
    parser.add_argument('--g4-file', required=True, help='G-quadruplex results CSV file')
    parser.add_argument('--gtf-file', required=True, help='GTF annotation file')
    parser.add_argument('--output-dir', default='results', help='Output directory')
    parser.add_argument('--promoter-upstream', type=int, default=1000, 
                       help='Promoter upstream region (bp)')
    parser.add_argument('--promoter-downstream', type=int, default=1000,
                       help='Promoter downstream region (bp)')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("=== Genomic Location Analysis ===")
    print(f"Z-DNA file: {args.zdna_file}")
    print(f"G4 file: {args.g4_file}")
    print(f"GTF file: {args.gtf_file}")
    print(f"Output directory: {args.output_dir}")
    
    # Load data
    print("\nLoading data...")
    genes_df = parse_gtf_file(args.gtf_file)
    
    print("Creating promoter regions...")
    promoters_df = create_promoter_regions(genes_df, args.promoter_upstream, args.promoter_downstream)
    
    print(f"Loaded {len(genes_df)} genes and created {len(promoters_df)} promoter regions")
    
    all_structures = []
    
    # Load Z-DNA results if available
    if args.zdna_file and Path(args.zdna_file).exists():
        print(f"📖 Loading Z-DNA results from {args.zdna_file}...")
        zdna_df = pd.read_csv(args.zdna_file)
        zdna_df['type'] = 'Z-DNA'
        all_structures.append(zdna_df)
        print(f"✅ Loaded {len(zdna_df)} Z-DNA structures")
    else:
        print("⚠️  No Z-DNA file provided or file not found, skipping Z-DNA analysis")
    
    # Load G4 results
    print(f"📖 Loading G-quadruplex results from {args.g4_file}...")
    g4_df = pd.read_csv(args.g4_file)
    g4_df['type'] = 'G4'
    all_structures.append(g4_df)
    print(f"✅ Loaded {len(g4_df)} G-quadruplex structures")
    
    # Combine all structures
    if all_structures:
        combined_df = pd.concat(all_structures, ignore_index=True)
        print(f"📊 Total structures: {len(combined_df)}")
        
        # Analyze genomic distribution
        chrom_dist = analyze_genomic_distribution(combined_df, genes_df, output_dir)
        
        # Find promoter overlaps
        promoter_overlaps = find_overlaps(
            combined_df, promoters_df, "structure"
        )
        
        print(f"🎯 Found {len(promoter_overlaps)} promoter overlaps")
        
        # Save results
        promoter_overlaps.to_csv(output_dir / 'promoter_overlaps.csv', index=False)
        
        # Create gene lists for STRING
        if len(promoter_overlaps) > 0:
            unique_genes = promoter_overlaps['gene_name'].unique()
            
            # Save gene list for STRING
            with open(output_dir / 'genes_for_string.txt', 'w') as f:
                for gene in unique_genes:
                    f.write(f"{gene}\n")
            
            print(f"💾 Saved {len(unique_genes)} unique genes to genes_for_string.txt")
            print("🔗 Upload this file to STRING-DB for functional enrichment analysis")
            
            # Create summary report
            summary = {
                'Total structures': len(combined_df),
                'Structures in promoters': len(promoter_overlaps),
                'Unique genes with structures': len(unique_genes),
                'Chromosomes analyzed': len(chrom_dist)
            }
            
            if args.zdna_file and Path(args.zdna_file).exists():
                zdna_promoter = promoter_overlaps[promoter_overlaps['feature_type'] == 'Z-DNA']
                g4_promoter = promoter_overlaps[promoter_overlaps['feature_type'] == 'G4']
                summary['Z-DNA in promoters'] = len(zdna_promoter)
                summary['G4 in promoters'] = len(g4_promoter)
            else:
                summary['G4 in promoters'] = len(promoter_overlaps)
            
            # Save summary
            with open(output_dir / 'analysis_summary.txt', 'w') as f:
                f.write("GENOMIC ANALYSIS SUMMARY\n")
                f.write("=" * 30 + "\n\n")
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")
            
            print("\n" + "=" * 50)
            print("📋 ANALYSIS SUMMARY:")
            for key, value in summary.items():
                print(f"   {key}: {value}")
            print("=" * 50)
            
        else:
            print("⚠️  No promoter overlaps found!")
    
    else:
        print("❌ No structure data loaded!")
        return 1
    
    print("✅ Analysis complete! Check results/ folder for outputs.")
    return 0

if __name__ == "__main__":
    sys.exit(main()) 