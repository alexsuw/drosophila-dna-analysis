#!/usr/bin/env python3
"""
Progress Monitor Script
Monitors the progress of Z-Hunt and G-quadruplex analysis
"""

import os
import time
import subprocess
import argparse

def check_process_status():
    """
    Check if Z-Hunt and G-quadruplex processes are running
    """
    try:
        # Check for zhunt process
        zhunt_result = subprocess.run(['pgrep', '-f', 'zhunt'], 
                                    capture_output=True, text=True)
        zhunt_running = bool(zhunt_result.stdout.strip())
        
        # Check for quadruplex search process
        python_result = subprocess.run(['pgrep', '-f', 'quadruplex_search'], 
                                     capture_output=True, text=True)
        python_running = bool(python_result.stdout.strip())
        
        return zhunt_running, python_running
    except:
        return False, False

def get_file_size(filepath):
    """
    Get file size in MB
    """
    try:
        size_bytes = os.path.getsize(filepath)
        return size_bytes / (1024 * 1024)  # Convert to MB
    except:
        return 0

def count_lines(filepath):
    """
    Count lines in a file
    """
    try:
        with open(filepath, 'r') as f:
            return sum(1 for line in f)
    except:
        return 0

def monitor_progress(interval=30):
    """
    Monitor progress of analysis
    """
    print("=== Analysis Progress Monitor ===")
    print("Monitoring Z-Hunt and G-quadruplex search...")
    print("Press Ctrl+C to stop monitoring\n")
    
    try:
        while True:
            zhunt_running, python_running = check_process_status()
            
            print(f"[{time.strftime('%H:%M:%S')}] Status:")
            
            # Z-Hunt status
            zdna_file = "data/results/z_dna_raw.txt"
            zdna_size = get_file_size(zdna_file)
            zdna_lines = count_lines(zdna_file)
            
            print(f"  Z-Hunt: {'🔄 Running' if zhunt_running else '✅ Completed'}")
            print(f"    Output file: {zdna_size:.1f} MB, {zdna_lines} lines")
            
            # G-quadruplex status
            g4_file = "results/quadruplex_results.csv"
            if os.path.exists(g4_file):
                g4_size = get_file_size(g4_file)
                g4_lines = count_lines(g4_file)
                print(f"  G-quadruplex: ✅ Completed")
                print(f"    Output file: {g4_size:.1f} MB, {g4_lines} lines")
            else:
                print(f"  G-quadruplex: {'🔄 Running' if python_running else '⏳ Waiting'}")
            
            # Check if both are done
            if not zhunt_running and not python_running and os.path.exists(g4_file):
                print("\n🎉 Both analyses completed!")
                print("\nNext steps:")
                print("1. Analyze Z-Hunt results")
                print("2. Analyze genomic locations")
                print("3. Perform functional enrichment")
                break
            
            print("-" * 50)
            time.sleep(interval)
            
    except KeyboardInterrupt:
        print("\nMonitoring stopped by user.")

def main():
    parser = argparse.ArgumentParser(description='Monitor analysis progress')
    parser.add_argument('--interval', type=int, default=30, 
                       help='Monitoring interval in seconds (default: 30)')
    
    args = parser.parse_args()
    monitor_progress(args.interval)

if __name__ == "__main__":
    main() 