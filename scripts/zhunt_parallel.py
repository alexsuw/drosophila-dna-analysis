#!/usr/bin/env python3
"""
Parallel Z-Hunt Analysis by Chromosome
"""

import subprocess
import time
import os
import sys
from pathlib import Path
import concurrent.futures
from collections import defaultdict

def split_genome_by_chromosome(fasta_file, output_dir):
    """Split genome FASTA by chromosomes"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"📄 Splitting {fasta_file} by chromosomes...")
    
    chromosomes = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous chromosome
                if current_chr:
                    chr_file = output_dir / f"{current_chr}.fa"
                    with open(chr_file, 'w') as cf:
                        cf.write(f">{current_chr}\n")
                        cf.write(''.join(current_seq) + '\n')
                    chromosomes[current_chr] = str(chr_file)
                
                # Start new chromosome
                current_chr = line[1:].split()[0]  # Get chromosome name
                current_seq = []
                print(f"   Found chromosome: {current_chr}")
            else:
                current_seq.append(line)
        
        # Save last chromosome
        if current_chr:
            chr_file = output_dir / f"{current_chr}.fa"
            with open(chr_file, 'w') as cf:
                cf.write(f">{current_chr}\n")
                cf.write(''.join(current_seq) + '\n')
            chromosomes[current_chr] = str(chr_file)
    
    print(f"✅ Split into {len(chromosomes)} chromosomes")
    return chromosomes

def run_zhunt_on_chromosome(chr_name, chr_file, output_dir, params=(12, 8, 12)):
    """Run Z-Hunt on a single chromosome"""
    output_file = Path(output_dir) / f"{chr_name}_zhunt.txt"
    
    start_time = time.time()
    cmd = ["./tools/zhunt/zhunt2", str(params[0]), str(params[1]), str(params[2]), chr_file]
    
    print(f"🔬 Starting {chr_name}...")
    
    try:
        with open(output_file, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True)
        
        elapsed = time.time() - start_time
        file_size = os.path.getsize(output_file) / (1024 * 1024)
        
        if result.returncode == 0:
            # Count lines
            with open(output_file, 'r') as f:
                line_count = sum(1 for line in f if line.strip())
            
            print(f"✅ {chr_name}: {elapsed:.1f}s, {file_size:.1f}MB, {line_count} regions")
            return {
                'chromosome': chr_name,
                'success': True,
                'time': elapsed,
                'size': file_size,
                'lines': line_count,
                'output_file': str(output_file)
            }
        else:
            print(f"❌ {chr_name}: Failed - {result.stderr}")
            return {
                'chromosome': chr_name,
                'success': False,
                'error': result.stderr,
                'output_file': str(output_file)
            }
            
    except Exception as e:
        print(f"❌ {chr_name}: Exception - {e}")
        return {
            'chromosome': chr_name,
            'success': False,
            'error': str(e),
            'output_file': str(output_file)
        }

def combine_results(results, output_file):
    """Combine chromosome results into single file"""
    print(f"🔗 Combining results into {output_file}...")
    
    total_lines = 0
    with open(output_file, 'w') as outf:
        for result in results:
            if result['success'] and os.path.exists(result['output_file']):
                with open(result['output_file'], 'r') as inf:
                    for line in inf:
                        if line.strip():
                            outf.write(line)
                            total_lines += 1
    
    file_size = os.path.getsize(output_file) / (1024 * 1024)
    print(f"✅ Combined: {total_lines:,} lines, {file_size:.1f} MB")
    return total_lines, file_size

def main():
    if len(sys.argv) != 3:
        print("Usage: python zhunt_parallel.py <genome.fa> <output.txt>")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(genome_file):
        print(f"❌ Genome file not found: {genome_file}")
        sys.exit(1)
    
    print("=" * 60)
    print("🧬 Parallel Z-Hunt Analysis by Chromosome")
    print("=" * 60)
    
    # Create working directory
    work_dir = Path("temp_chromosomes")
    work_dir.mkdir(exist_ok=True)
    
    # Split genome
    start_time = time.time()
    chromosomes = split_genome_by_chromosome(genome_file, work_dir)
    
    # Filter out small chromosomes (< 1MB) to speed up
    large_chromosomes = {}
    for chr_name, chr_file in chromosomes.items():
        size = os.path.getsize(chr_file) / (1024 * 1024)
        if size > 1.0:  # Only process chromosomes > 1MB
            large_chromosomes[chr_name] = chr_file
        else:
            print(f"⏭️  Skipping small chromosome {chr_name} ({size:.1f}MB)")
    
    print(f"\n🚀 Processing {len(large_chromosomes)} large chromosomes in parallel...")
    
    # Run Z-Hunt in parallel (max 4 processes to avoid overwhelming system)
    results = []
    max_workers = min(4, len(large_chromosomes))
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_chr = {
            executor.submit(run_zhunt_on_chromosome, chr_name, chr_file, work_dir): chr_name
            for chr_name, chr_file in large_chromosomes.items()
        }
        
        # Collect results
        for future in concurrent.futures.as_completed(future_to_chr):
            result = future.result()
            results.append(result)
    
    # Combine successful results
    successful_results = [r for r in results if r['success']]
    if successful_results:
        total_lines, final_size = combine_results(successful_results, output_file)
        
        elapsed = time.time() - start_time
        
        print(f"\n🎉 Analysis complete!")
        print(f"   ⏱️  Total time: {elapsed:.1f}s")
        print(f"   📄 Output: {output_file}")
        print(f"   📊 Results: {total_lines:,} Z-DNA regions")
        print(f"   💾 File size: {final_size:.1f} MB")
        print(f"   ✅ Processed: {len(successful_results)}/{len(large_chromosomes)} chromosomes")
        
        # Cleanup
        print(f"\n🧹 Cleaning up temporary files...")
        import shutil
        shutil.rmtree(work_dir)
        
    else:
        print(f"\n❌ No successful results!")
        sys.exit(1)

if __name__ == "__main__":
    main() 