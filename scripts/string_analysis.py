#!/usr/bin/env python3
"""
STRING Database Functional Enrichment Analysis
"""

import requests
import pandas as pd
import json
import time
import sys
from pathlib import Path

def read_gene_list(gene_file):
    """Read gene list from file"""
    genes = []
    with open(gene_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.append(gene)
    return genes

def convert_flybase_to_string(flybase_genes, species_id=7227):
    """Convert FlyBase IDs to STRING identifiers"""
    print(f"🔗 Converting {len(flybase_genes)} FlyBase IDs to STRING format...")
    
    # STRING API endpoint for identifier mapping
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv"
    method = "get_string_ids"
    
    # Prepare parameters
    params = {
        "identifiers": "\r".join(flybase_genes),  # Join genes with carriage return
        "species": species_id,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "bioinfo_homework"
    }
    
    try:
        response = requests.post(f"{string_api_url}/{output_format}/{method}", data=params)
        response.raise_for_status()
        
        # Parse response
        lines = response.text.strip().split('\n')
        header = lines[0].split('\t')
        
        string_ids = []
        mapping = {}
        
        for line in lines[1:]:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    query_id = parts[0]
                    string_id = parts[1]
                    if string_id != 'Error':
                        string_ids.append(string_id)
                        mapping[query_id] = string_id
        
        print(f"✅ Mapped {len(string_ids)} genes to STRING IDs")
        return string_ids, mapping
        
    except Exception as e:
        print(f"❌ Error mapping genes: {e}")
        return [], {}

def get_functional_enrichment(string_ids, species_id=7227):
    """Get functional enrichment from STRING"""
    if not string_ids:
        print("❌ No STRING IDs to analyze")
        return None
    
    print(f"🧬 Getting functional enrichment for {len(string_ids)} genes...")
    
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv"
    method = "enrichment"
    
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species_id,
        "caller_identity": "bioinfo_homework"
    }
    
    try:
        response = requests.post(f"{string_api_url}/{output_format}/{method}", data=params)
        response.raise_for_status()
        
        # Parse enrichment results
        lines = response.text.strip().split('\n')
        if len(lines) < 2:
            print("⚠️  No enrichment results returned")
            return None
            
        header = lines[0].split('\t')
        
        enrichment_data = []
        for line in lines[1:]:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= len(header):
                    enrichment_data.append(dict(zip(header, parts)))
        
        df = pd.DataFrame(enrichment_data)
        
        # Convert numeric columns
        numeric_cols = ['number_of_genes', 'number_of_genes_in_background', 'pvalue', 'fdr']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        print(f"✅ Found {len(df)} enrichment terms")
        return df
        
    except Exception as e:
        print(f"❌ Error getting enrichment: {e}")
        return None

def get_protein_interactions(string_ids, species_id=7227):
    """Get protein-protein interactions"""
    if not string_ids:
        return None
        
    print(f"🕸️  Getting protein interactions...")
    
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv"
    method = "network"
    
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species_id,
        "caller_identity": "bioinfo_homework"
    }
    
    try:
        response = requests.post(f"{string_api_url}/{output_format}/{method}", data=params)
        response.raise_for_status()
        
        lines = response.text.strip().split('\n')
        if len(lines) < 2:
            return None
            
        header = lines[0].split('\t')
        
        network_data = []
        for line in lines[1:]:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= len(header):
                    network_data.append(dict(zip(header, parts)))
        
        df = pd.DataFrame(network_data)
        
        # Convert score column
        if 'score' in df.columns:
            df['score'] = pd.to_numeric(df['score'], errors='coerce')
        
        print(f"✅ Found {len(df)} protein interactions")
        return df
        
    except Exception as e:
        print(f"❌ Error getting interactions: {e}")
        return None

def save_results(enrichment_df, network_df, mapping, output_dir):
    """Save all results"""
    output_dir = Path(output_dir)
    
    print(f"💾 Saving results to {output_dir}/...")
    
    # Save enrichment results
    if enrichment_df is not None and not enrichment_df.empty:
        # Filter significant results (FDR < 0.05)
        significant = enrichment_df[enrichment_df['fdr'] < 0.05].copy()
        significant = significant.sort_values('fdr')
        
        enrichment_df.to_csv(output_dir / 'string_enrichment_all.csv', index=False)
        significant.to_csv(output_dir / 'string_enrichment_significant.csv', index=False)
        
        print(f"   📊 Enrichment: {len(enrichment_df)} total, {len(significant)} significant (FDR<0.05)")
        
        # Top categories
        if len(significant) > 0:
            print("   🏆 Top significant categories:")
            for _, row in significant.head(5).iterrows():
                print(f"      {row['category']}: {row['description']} (FDR={row['fdr']:.2e})")
    
    # Save network
    if network_df is not None and not network_df.empty:
        network_df.to_csv(output_dir / 'string_network.csv', index=False)
        print(f"   🕸️  Network: {len(network_df)} interactions")
    
    # Save mapping
    with open(output_dir / 'string_gene_mapping.json', 'w') as f:
        json.dump(mapping, f, indent=2)
    
    # Create summary
    summary = {
        'analysis_type': 'STRING Database Functional Enrichment',
        'species': 'Drosophila melanogaster (7227)',
        'input_genes': len(mapping),
        'mapped_genes': len([k for k in mapping.values() if k]),
        'enrichment_terms': len(enrichment_df) if enrichment_df is not None else 0,
        'significant_terms': len(enrichment_df[enrichment_df['fdr'] < 0.05]) if enrichment_df is not None else 0,
        'protein_interactions': len(network_df) if network_df is not None else 0
    }
    
    with open(output_dir / 'string_analysis_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✅ STRING analysis complete!")
    return summary

def main():
    if len(sys.argv) != 3:
        print("Usage: python string_analysis.py <gene_list.txt> <output_dir>")
        sys.exit(1)
    
    gene_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not Path(gene_file).exists():
        print(f"❌ Gene file not found: {gene_file}")
        sys.exit(1)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    print("=" * 60)
    print("🧬 STRING Database Functional Enrichment Analysis")
    print("=" * 60)
    print(f"Gene list: {gene_file}")
    print(f"Output: {output_dir}")
    print(f"Species: Drosophila melanogaster")
    
    # Read genes
    genes = read_gene_list(gene_file)
    print(f"📖 Loaded {len(genes)} genes")
    
    # Convert to STRING IDs
    string_ids, mapping = convert_flybase_to_string(genes)
    
    if not string_ids:
        print("❌ No genes could be mapped to STRING")
        sys.exit(1)
    
    # Get functional enrichment
    print("\n📊 Running functional enrichment analysis...")
    enrichment_df = get_functional_enrichment(string_ids)
    
    # Get protein interactions
    print("\n🕸️  Getting protein interaction network...")
    network_df = get_protein_interactions(string_ids)
    
    # Save results
    print("\n💾 Saving results...")
    summary = save_results(enrichment_df, network_df, mapping, output_dir)
    
    print(f"\n🎉 Analysis complete!")
    print(f"   📊 {summary['significant_terms']} significant GO terms found")
    print(f"   🕸️  {summary['protein_interactions']} protein interactions")
    print(f"\n📁 Check {output_dir}/ for detailed results")

if __name__ == "__main__":
    main() 