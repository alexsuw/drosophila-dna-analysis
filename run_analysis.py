#!/usr/bin/env python3
"""
Master Analysis Pipeline Script
Runs the complete Z-DNA and G-quadruplex analysis pipeline
"""

import os
import sys
import subprocess
import argparse
import time

def run_command(cmd, description, background=False):
    """
    Run a command with error handling
    """
    print(f"\n{'='*50}")
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    print(f"{'='*50}")
    
    if background:
        process = subprocess.Popen(cmd, shell=True)
        print(f"Started in background (PID: {process.pid})")
        return process
    else:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR: {description} failed!")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return False
        else:
            print(f"SUCCESS: {description} completed!")
            if result.stdout:
                print(f"Output: {result.stdout}")
            return True

def check_file_exists(filepath, description):
    """
    Check if a required file exists
    """
    if not os.path.exists(filepath):
        print(f"ERROR: {description} not found at {filepath}")
        return False
    print(f"✓ Found {description}: {filepath}")
    return True

def wait_for_completion(processes, descriptions):
    """
    Wait for background processes to complete
    """
    print(f"\nWaiting for {len(processes)} processes to complete...")
    for i, (process, desc) in enumerate(zip(processes, descriptions)):
        print(f"Waiting for: {desc}")
        process.wait()
        if process.returncode == 0:
            print(f"✓ {desc} completed successfully")
        else:
            print(f"✗ {desc} failed with return code {process.returncode}")

def main():
    parser = argparse.ArgumentParser(description='Run complete Z-DNA and G-quadruplex analysis')
    parser.add_argument('--skip-download', action='store_true', 
                       help='Skip genome download (if already downloaded)')
    parser.add_argument('--skip-zhunt', action='store_true',
                       help='Skip Z-Hunt analysis (if already completed)')
    parser.add_argument('--skip-g4', action='store_true',
                       help='Skip G-quadruplex search (if already completed)')
    parser.add_argument('--min-zscore', type=float, default=300,
                       help='Minimum Z-score threshold for filtering')
    parser.add_argument('--max-zscore', type=float, default=400,
                       help='Maximum Z-score threshold for filtering')
    parser.add_argument('--g4-min-score', type=float, default=60,
                       help='Minimum G-quadruplex score threshold')
    
    args = parser.parse_args()
    
    print("🧬 Z-DNA and G-Quadruplex Analysis Pipeline")
    print("=" * 60)
    
    # Check prerequisites
    print("\n1. Checking prerequisites...")
    
    if not args.skip_download:
        # Check if genome files exist
        genome_file = "data/genome/dm6.fa"
        gtf_file = "data/annotation/dm6.ensGene.gtf"
        
        if not check_file_exists(genome_file, "Genome FASTA"):
            print("Please download the genome first or use --skip-download")
            return 1
        if not check_file_exists(gtf_file, "Gene annotation GTF"):
            print("Please download the annotation first or use --skip-download")
            return 1
    
    # Check Z-Hunt executable
    zhunt_exe = "tools/zhunt/zhunt2"
    if not check_file_exists(zhunt_exe, "Z-Hunt executable"):
        print("Please compile Z-Hunt first")
        return 1
    
    # Create output directories
    os.makedirs("results", exist_ok=True)
    os.makedirs("data/results", exist_ok=True)
    
    background_processes = []
    process_descriptions = []
    
    # Step 2: Run Z-Hunt analysis
    if not args.skip_zhunt:
        print("\n2. Running Z-Hunt analysis...")
        zhunt_cmd = f"./tools/zhunt/zhunt2 12 8 12 data/genome/dm6.fa > data/results/z_dna_raw.txt"
        process = run_command(zhunt_cmd, "Z-Hunt analysis", background=True)
        if process:
            background_processes.append(process)
            process_descriptions.append("Z-Hunt analysis")
    else:
        print("\n2. Skipping Z-Hunt analysis (--skip-zhunt specified)")
    
    # Step 3: Run G-quadruplex search
    if not args.skip_g4:
        print("\n3. Running G-quadruplex search...")
        g4_cmd = f"python scripts/quadruplex_search.py --input data/genome/dm6.fa --output-dir results --min-score {args.g4_min_score}"
        process = run_command(g4_cmd, "G-quadruplex search", background=True)
        if process:
            background_processes.append(process)
            process_descriptions.append("G-quadruplex search")
    else:
        print("\n3. Skipping G-quadruplex search (--skip-g4 specified)")
    
    # Wait for background processes
    if background_processes:
        wait_for_completion(background_processes, process_descriptions)
    
    # Step 4: Analyze Z-Hunt results
    print("\n4. Analyzing Z-Hunt results...")
    zhunt_analysis_cmd = f"python scripts/zhunt_analysis.py --input data/results/z_dna_raw.txt --output-dir results --min-zscore {args.min_zscore} --max-zscore {args.max_zscore}"
    if not run_command(zhunt_analysis_cmd, "Z-Hunt results analysis"):
        return 1
    
    # Step 5: Genomic location analysis
    print("\n5. Running genomic location analysis...")
    location_cmd = f"python scripts/analysis.py --zdna-file results/zdna_filtered.csv --g4-file results/quadruplex_results.csv --gtf-file data/annotation/dm6.ensGene.gtf --output-dir results"
    if not run_command(location_cmd, "Genomic location analysis"):
        return 1
    
    # Step 6: Summary
    print("\n" + "="*60)
    print("🎉 ANALYSIS PIPELINE COMPLETED!")
    print("="*60)
    
    print("\nGenerated files:")
    result_files = [
        "results/zdna_filtered.csv",
        "results/zdna_filtered.bed", 
        "results/quadruplex_results.csv",
        "results/quadruplex_results.bed",
        "results/zdna_analysis.png",
        "results/quadruplex_analysis.png",
        "results/genomic_analysis.png",
        "results/zdna_genes.txt",
        "results/zdna_promoter_genes.txt",
        "results/g4_genes.txt",
        "results/g4_promoter_genes.txt"
    ]
    
    for file in result_files:
        if os.path.exists(file):
            size = os.path.getsize(file) / 1024  # KB
            print(f"  ✓ {file} ({size:.1f} KB)")
        else:
            print(f"  ✗ {file} (not found)")
    
    print("\nNext steps for functional enrichment analysis:")
    print("1. Go to https://string-db.org/")
    print("2. Upload gene lists from results/ directory:")
    print("   - zdna_promoter_genes.txt")
    print("   - g4_promoter_genes.txt")
    print("3. Select 'Drosophila melanogaster' as organism")
    print("4. Analyze functional enrichment and pathways")
    print("5. Export results and interpret biological significance")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())