#!/usr/bin/env python3
"""
Z-Hunt with Progress Monitoring
"""

import subprocess
import time
import os
import sys
from pathlib import Path

def get_file_size(filepath):
    """Get file size in MB"""
    try:
        return os.path.getsize(filepath) / (1024 * 1024)
    except:
        return 0

def get_genome_size(fasta_file):
    """Estimate genome size for progress calculation"""
    print("📏 Calculating genome size...")
    size = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                size += len(line.strip())
    return size

def run_zhunt_with_progress(genome_file, output_file, min_size=12, window=8, max_size=12):
    """Run Z-Hunt with progress monitoring"""
    
    # Create output directory
    output_dir = Path(output_file).parent
    output_dir.mkdir(exist_ok=True)
    
    print(f"🔬 Starting Z-Hunt analysis...")
    print(f"   Genome: {genome_file}")
    print(f"   Output: {output_file}")
    print(f"   Parameters: {min_size} {window} {max_size}")
    
    # Get genome size for progress estimation
    genome_size = get_genome_size(genome_file)
    print(f"   Genome size: {genome_size:,} bp")
    
    # Start Z-Hunt process
    cmd = ["./tools/zhunt/zhunt2", str(min_size), str(window), str(max_size), genome_file]
    
    print(f"\n🚀 Running: {' '.join(cmd)}")
    print("📊 Progress monitoring:")
    
    start_time = time.time()
    
    # Run process and redirect output
    with open(output_file, 'w') as outf:
        process = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.PIPE, text=True)
        
        # Monitor progress
        while process.poll() is None:
            # Check output file size
            file_size = get_file_size(output_file)
            elapsed = time.time() - start_time
            
            # Estimate progress (very rough)
            # Z-Hunt typically produces ~1-10 MB per million bp
            estimated_progress = min(95, (file_size / (genome_size / 1000000)) * 10)
            
            print(f"\r   ⏱️  {elapsed:.0f}s | 📄 {file_size:.1f} MB | 📈 {estimated_progress:.1f}% | Status: Running...", 
                  end="", flush=True)
            
            time.sleep(5)  # Update every 5 seconds
        
        # Process finished
        return_code = process.returncode
        stderr = process.stderr.read()
        
        elapsed = time.time() - start_time
        final_size = get_file_size(output_file)
        
        print(f"\r   ⏱️  {elapsed:.0f}s | 📄 {final_size:.1f} MB | ✅ Complete!                    ")
        
        if return_code == 0:
            print(f"🎉 Z-Hunt completed successfully!")
            
            # Count lines in output
            try:
                with open(output_file, 'r') as f:
                    line_count = sum(1 for line in f if line.strip())
                print(f"   📊 Found {line_count:,} potential Z-DNA regions")
            except:
                print("   ⚠️  Could not count results")
                
        else:
            print(f"❌ Z-Hunt failed with return code: {return_code}")
            if stderr:
                print(f"   Error: {stderr}")
                
        return return_code, final_size

def main():
    if len(sys.argv) != 3:
        print("Usage: python zhunt_with_progress.py <genome.fa> <output.txt>")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(genome_file):
        print(f"❌ Genome file not found: {genome_file}")
        sys.exit(1)
    
    print("=" * 60)
    print("🧬 Z-Hunt Analysis with Progress Monitoring")
    print("=" * 60)
    
    return_code, final_size = run_zhunt_with_progress(genome_file, output_file)
    
    if return_code == 0:
        print(f"\n✅ Analysis complete! Output saved to: {output_file}")
        print(f"   Final file size: {final_size:.1f} MB")
        print(f"\n🔗 Next steps:")
        print(f"   1. Process results with: python scripts/zhunt_analysis.py")
        print(f"   2. Run genomic analysis: python scripts/analysis.py")
    else:
        print(f"\n❌ Analysis failed!")
        sys.exit(1)

if __name__ == "__main__":
    main() 