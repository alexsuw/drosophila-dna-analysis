#!/usr/bin/env python3
"""
Comprehensive Visualization Script for G-quadruplex and STRING Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pathlib import Path
import sys

# Set style for beautiful plots
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

def load_data(results_dir):
    """Load all analysis results"""
    results_dir = Path(results_dir)
    
    data = {}
    
    # Load G-quadruplex results
    g4_file = results_dir / 'quadruplex_results.csv'
    if g4_file.exists():
        data['g4'] = pd.read_csv(g4_file)
        print(f"✅ Loaded {len(data['g4'])} G-quadruplex structures")
    
    # Load promoter overlaps
    promoter_file = results_dir / 'promoter_overlaps.csv'
    if promoter_file.exists():
        data['promoters'] = pd.read_csv(promoter_file)
        print(f"✅ Loaded {len(data['promoters'])} promoter overlaps")
    
    # Load STRING enrichment
    string_file = results_dir / 'string_enrichment_significant.csv'
    if string_file.exists():
        data['string'] = pd.read_csv(string_file)
        print(f"✅ Loaded {len(data['string'])} STRING enrichment terms")
    
    # Load STRING network
    network_file = results_dir / 'string_network.csv'
    if network_file.exists():
        data['network'] = pd.read_csv(network_file)
        print(f"✅ Loaded {len(data['network'])} protein interactions")
    
    return data

def create_g4_overview(data, output_dir):
    """Create G-quadruplex overview plots"""
    g4_df = data['g4']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('G-Quadruplex Analysis Overview', fontsize=16, fontweight='bold')
    
    # 1. Chromosome distribution
    chrom_counts = g4_df['chromosome'].value_counts().head(10)
    axes[0, 0].bar(range(len(chrom_counts)), chrom_counts.values, 
                   color=sns.color_palette("viridis", len(chrom_counts)))
    axes[0, 0].set_xticks(range(len(chrom_counts)))
    axes[0, 0].set_xticklabels(chrom_counts.index, rotation=45)
    axes[0, 0].set_title('G4 Distribution by Chromosome')
    axes[0, 0].set_ylabel('Number of G4 structures')
    
    # Add value labels
    for i, v in enumerate(chrom_counts.values):
        axes[0, 0].text(i, v + 20, str(v), ha='center', va='bottom')
    
    # 2. Length distribution
    axes[0, 1].hist(g4_df['length'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 1].axvline(g4_df['length'].mean(), color='red', linestyle='--', 
                       label=f'Mean: {g4_df["length"].mean():.1f} bp')
    axes[0, 1].set_title('G4 Length Distribution')
    axes[0, 1].set_xlabel('Length (bp)')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].legend()
    
    # 3. Score distribution
    axes[0, 2].hist(g4_df['score'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    axes[0, 2].axvline(g4_df['score'].mean(), color='red', linestyle='--',
                       label=f'Mean: {g4_df["score"].mean():.1f}')
    axes[0, 2].set_title('G4 Score Distribution')
    axes[0, 2].set_xlabel('G4 Score')
    axes[0, 2].set_ylabel('Frequency')
    axes[0, 2].legend()
    
    # 4. GC content vs Score
    axes[1, 0].scatter(g4_df['gc_content'], g4_df['score'], alpha=0.6, color='purple')
    axes[1, 0].set_title('GC Content vs G4 Score')
    axes[1, 0].set_xlabel('GC Content')
    axes[1, 0].set_ylabel('G4 Score')
    
    # Add correlation
    corr = g4_df['gc_content'].corr(g4_df['score'])
    axes[1, 0].text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                    transform=axes[1, 0].transAxes, bbox=dict(boxstyle="round", facecolor='white'))
    
    # 5. G-run length distribution
    g_run_counts = g4_df['g_run_length'].value_counts().sort_index()
    axes[1, 1].bar(g_run_counts.index, g_run_counts.values, color='orange', alpha=0.7)
    axes[1, 1].set_title('G-run Length Distribution')
    axes[1, 1].set_xlabel('G-run Length')
    axes[1, 1].set_ylabel('Number of G4 structures')
    
    # 6. Summary statistics
    axes[1, 2].axis('off')
    stats_text = f"""
    G-Quadruplex Summary Statistics
    
    Total G4 structures: {len(g4_df):,}
    
    Length:
    • Mean: {g4_df['length'].mean():.1f} bp
    • Median: {g4_df['length'].median():.1f} bp
    • Range: {g4_df['length'].min()}-{g4_df['length'].max()} bp
    
    Score:
    • Mean: {g4_df['score'].mean():.1f}
    • Median: {g4_df['score'].median():.1f}
    • Max: {g4_df['score'].max():.1f}
    
    GC Content:
    • Mean: {g4_df['gc_content'].mean():.3f}
    • Range: {g4_df['gc_content'].min():.3f}-{g4_df['gc_content'].max():.3f}
    
    Chromosomes: {g4_df['chromosome'].nunique()}
    """
    axes[1, 2].text(0.1, 0.9, stats_text, transform=axes[1, 2].transAxes, 
                     fontsize=11, verticalalignment='top',
                     bbox=dict(boxstyle="round,pad=0.5", facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'g4_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ G4 overview plots created")

def create_promoter_analysis(data, output_dir):
    """Create promoter analysis plots"""
    if 'promoters' not in data or data['promoters'].empty:
        print("⚠️  No promoter data available")
        return
    
    promoters_df = data['promoters']
    g4_df = data['g4']
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('G-Quadruplex Promoter Analysis', fontsize=16, fontweight='bold')
    
    # 1. G4 score distribution in promoters
    axes[0, 0].hist(promoters_df['score'], bins=30, alpha=0.7, 
                    color='lightcoral', edgecolor='black')
    axes[0, 0].set_title('G4 Score Distribution in Promoters')
    axes[0, 0].set_xlabel('G4 Score')
    axes[0, 0].set_ylabel('Number of G4 structures')
    axes[0, 0].axvline(promoters_df['score'].mean(), color='red', linestyle='--',
                       label=f'Mean: {promoters_df["score"].mean():.1f}')
    axes[0, 0].legend()
    
    # 2. Promoter vs non-promoter G4s
    total_g4 = len(g4_df)
    promoter_g4 = len(promoters_df)
    non_promoter_g4 = total_g4 - promoter_g4
    
    labels = ['In Promoters', 'Outside Promoters']
    sizes = [promoter_g4, non_promoter_g4]
    colors = ['lightgreen', 'lightgray']
    
    wedges, texts, autotexts = axes[0, 1].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                                              startangle=90, explode=(0.05, 0))
    axes[0, 1].set_title('G4 Distribution: Promoters vs Genome')
    
    # 3. Chromosome distribution of promoter G4s
    prom_chrom = promoters_df['feature_chr'].value_counts().head(8)
    axes[1, 0].bar(range(len(prom_chrom)), prom_chrom.values, 
                   color=sns.color_palette("Set2", len(prom_chrom)))
    axes[1, 0].set_xticks(range(len(prom_chrom)))
    axes[1, 0].set_xticklabels(prom_chrom.index, rotation=45)
    axes[1, 0].set_title('Promoter G4s by Chromosome')
    axes[1, 0].set_ylabel('Number of G4 structures')
    
    # 4. Feature type distribution
    if 'feature_type' in promoters_df.columns:
        type_counts = promoters_df['feature_type'].value_counts()
        axes[1, 1].bar(range(len(type_counts)), type_counts.values, 
                       color=sns.color_palette("Set1", len(type_counts)))
        axes[1, 1].set_xticks(range(len(type_counts)))
        axes[1, 1].set_xticklabels(type_counts.index, rotation=45)
        axes[1, 1].set_title('G4s by Feature Type')
        axes[1, 1].set_xlabel('Feature Type')
        axes[1, 1].set_ylabel('Number of G4 structures')
        
        # Add percentages
        for i, count in enumerate(type_counts.values):
            pct = count / type_counts.sum() * 100
            axes[1, 1].text(i, count + 10, f'{pct:.1f}%', ha='center', va='bottom')
    else:
        axes[1, 1].axis('off')
        axes[1, 1].text(0.5, 0.5, 'Feature type data not available', 
                         ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'promoter_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Promoter analysis plots created")

def create_string_enrichment_plots(data, output_dir):
    """Create STRING enrichment visualization"""
    if 'string' not in data:
        print("⚠️  No STRING data available")
        return
    
    string_df = data['string']
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('STRING Database Functional Enrichment Analysis', fontsize=16, fontweight='bold')
    
    # 1. Top enriched categories
    top_terms = string_df.nsmallest(15, 'fdr')
    
    y_pos = np.arange(len(top_terms))
    bars = axes[0, 0].barh(y_pos, -np.log10(top_terms['fdr']), 
                           color=sns.color_palette("viridis", len(top_terms)))
    axes[0, 0].set_yticks(y_pos)
    axes[0, 0].set_yticklabels([desc[:50] + '...' if len(desc) > 50 else desc 
                                for desc in top_terms['description']], fontsize=9)
    axes[0, 0].set_xlabel('-log10(FDR)')
    axes[0, 0].set_title('Top 15 Enriched Terms')
    axes[0, 0].invert_yaxis()
    
    # 2. Categories distribution
    category_counts = string_df['category'].value_counts()
    axes[0, 1].pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%',
                   colors=sns.color_palette("Set3", len(category_counts)))
    axes[0, 1].set_title('Enrichment by Category Type')
    
    # 3. Gene count vs significance
    axes[1, 0].scatter(string_df['number_of_genes'], -np.log10(string_df['fdr']), 
                       alpha=0.6, s=50, color='purple')
    axes[1, 0].set_xlabel('Number of Genes in Term')
    axes[1, 0].set_ylabel('-log10(FDR)')
    axes[1, 0].set_title('Gene Count vs Significance')
    
    # Add significance threshold line
    axes[1, 0].axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7,
                       label='FDR = 0.05')
    axes[1, 0].legend()
    
    # 4. Top biological processes
    bio_processes = string_df[string_df['category'] == 'Process'].nsmallest(10, 'fdr')
    if not bio_processes.empty:
        y_pos = np.arange(len(bio_processes))
        axes[1, 1].barh(y_pos, -np.log10(bio_processes['fdr']), color='lightgreen')
        axes[1, 1].set_yticks(y_pos)
        axes[1, 1].set_yticklabels([desc[:40] + '...' if len(desc) > 40 else desc 
                                    for desc in bio_processes['description']], fontsize=8)
        axes[1, 1].set_xlabel('-log10(FDR)')
        axes[1, 1].set_title('Top Biological Processes')
        axes[1, 1].invert_yaxis()
    else:
        axes[1, 1].text(0.5, 0.5, 'No biological processes found', 
                         ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'string_enrichment_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ STRING enrichment plots created")

def create_summary_dashboard(data, output_dir):
    """Create a comprehensive summary dashboard"""
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    fig.suptitle('G-Quadruplex & Functional Enrichment Analysis Dashboard', 
                 fontsize=18, fontweight='bold')
    
    # Main statistics panel
    ax_stats = fig.add_subplot(gs[0, :2])
    ax_stats.axis('off')
    
    g4_df = data.get('g4', pd.DataFrame())
    promoters_df = data.get('promoters', pd.DataFrame())
    string_df = data.get('string', pd.DataFrame())
    
    stats_text = f"""
    📊 ANALYSIS SUMMARY
    
    🧬 G-Quadruplex Structures:
    • Total identified: {len(g4_df):,}
    • Mean length: {g4_df['length'].mean():.1f} bp
    • Mean score: {g4_df['score'].mean():.1f}
    • Chromosomes covered: {g4_df['chromosome'].nunique() if not g4_df.empty else 0}
    
    🎯 Promoter Analysis:
    • G4s in promoters: {len(promoters_df):,}
    • Unique genes affected: {promoters_df['gene_name'].nunique() if not promoters_df.empty else 0}
    • Promoter enrichment: {len(promoters_df)/len(g4_df)*100:.1f}% of all G4s
    
    🔬 Functional Enrichment:
    • Significant GO terms: {len(string_df):,}
    • Most significant: {string_df.iloc[0]['description'][:50] + '...' if not string_df.empty else 'N/A'}
    • FDR range: {string_df['fdr'].min():.2e} - {string_df['fdr'].max():.2e}
    """
    
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, fontsize=12,
                  verticalalignment='top', bbox=dict(boxstyle="round,pad=0.5", 
                  facecolor='lightblue', alpha=0.8))
    
    # G4 chromosome distribution (compact)
    ax_chrom = fig.add_subplot(gs[0, 2:])
    if not g4_df.empty:
        chrom_counts = g4_df['chromosome'].value_counts().head(8)
        bars = ax_chrom.bar(range(len(chrom_counts)), chrom_counts.values, 
                           color=sns.color_palette("viridis", len(chrom_counts)))
        ax_chrom.set_xticks(range(len(chrom_counts)))
        ax_chrom.set_xticklabels(chrom_counts.index, rotation=45)
        ax_chrom.set_title('G4 Distribution by Chromosome')
        ax_chrom.set_ylabel('Count')
        
        # Add value labels
        for i, v in enumerate(chrom_counts.values):
            ax_chrom.text(i, v + 10, str(v), ha='center', va='bottom', fontsize=9)
    
    # Length and score distributions
    ax_length = fig.add_subplot(gs[1, 0])
    ax_score = fig.add_subplot(gs[1, 1])
    
    if not g4_df.empty:
        ax_length.hist(g4_df['length'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax_length.set_title('G4 Length Distribution')
        ax_length.set_xlabel('Length (bp)')
        ax_length.set_ylabel('Frequency')
        
        ax_score.hist(g4_df['score'], bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
        ax_score.set_title('G4 Score Distribution')
        ax_score.set_xlabel('Score')
        ax_score.set_ylabel('Frequency')
    
    # Promoter pie chart
    ax_pie = fig.add_subplot(gs[1, 2])
    if not g4_df.empty and not promoters_df.empty:
        total_g4 = len(g4_df)
        promoter_g4 = len(promoters_df)
        non_promoter_g4 = total_g4 - promoter_g4
        
        sizes = [promoter_g4, non_promoter_g4]
        labels = ['In Promoters', 'Outside Promoters']
        colors = ['lightcoral', 'lightgray']
        
        ax_pie.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax_pie.set_title('G4 Promoter Distribution')
    
    # Top STRING terms
    ax_string = fig.add_subplot(gs[1, 3])
    if not string_df.empty:
        top_5 = string_df.nsmallest(5, 'fdr')
        y_pos = np.arange(len(top_5))
        ax_string.barh(y_pos, -np.log10(top_5['fdr']), color='orange')
        ax_string.set_yticks(y_pos)
        ax_string.set_yticklabels([desc[:25] + '...' if len(desc) > 25 else desc 
                                   for desc in top_5['description']], fontsize=8)
        ax_string.set_xlabel('-log10(FDR)')
        ax_string.set_title('Top 5 GO Terms')
        ax_string.invert_yaxis()
    
    # Bottom row: detailed enrichment
    ax_enrich = fig.add_subplot(gs[2, :])
    if not string_df.empty:
        # Category breakdown
        category_counts = string_df['category'].value_counts()
        x_pos = np.arange(len(category_counts))
        bars = ax_enrich.bar(x_pos, category_counts.values, 
                            color=sns.color_palette("Set2", len(category_counts)))
        ax_enrich.set_xticks(x_pos)
        ax_enrich.set_xticklabels(category_counts.index, rotation=45)
        ax_enrich.set_title('Functional Categories Distribution')
        ax_enrich.set_ylabel('Number of Terms')
        
        # Add value labels
        for i, v in enumerate(category_counts.values):
            ax_enrich.text(i, v + 1, str(v), ha='center', va='bottom')
    
    plt.savefig(output_dir / 'comprehensive_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ Comprehensive dashboard created")

def main():
    if len(sys.argv) != 2:
        print("Usage: python create_visualizations.py <results_dir>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_dir = Path(results_dir)
    
    print("=" * 60)
    print("🎨 Creating Comprehensive Visualizations")
    print("=" * 60)
    
    # Load all data
    data = load_data(results_dir)
    
    if not data:
        print("❌ No data found to visualize")
        sys.exit(1)
    
    # Create visualizations
    print("\n🎨 Creating visualizations...")
    
    create_g4_overview(data, output_dir)
    create_promoter_analysis(data, output_dir)
    create_string_enrichment_plots(data, output_dir)
    create_summary_dashboard(data, output_dir)
    
    print(f"\n🎉 All visualizations created!")
    print(f"📁 Check {output_dir}/ for the following files:")
    print("   • g4_comprehensive_analysis.png")
    print("   • promoter_analysis.png") 
    print("   • string_enrichment_analysis.png")
    print("   • comprehensive_dashboard.png")

if __name__ == "__main__":
    main() 