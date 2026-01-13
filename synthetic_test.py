#!/usr/bin/env python3
# synthetic_test.py
# -*- coding: utf-8 -*-


import os
import re
import time
import subprocess


TEST_PARAMS = [
    "--beam_width", "20",      
    "--min_contig_len", "300"  
]

def get_baseline(error_rate_percent: float) -> float:
    """Return baseline score for given error rate.
    """
    if error_rate_percent <= 1.5:
        return 0.59
    elif error_rate_percent <= 3.5:
        return 0.21
    else:
        return 0.08


def get_process_memory(pid):
    try:
        with open(f"/proc/{pid}/status", "r") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return int(parts[1]) / 1024.0
    except (FileNotFoundError, ProcessLookupError):
        return 0.0

    return 0.0

def run_assembly_and_eval(input_file):
    output_file = input_file.replace(".fasta", "_contigs.fasta")
    cmd_asm = ["python3", "assembly.py", input_file, output_file] + TEST_PARAMS 
    start_time = time.time()
    max_mem = 0.0
    process = subprocess.Popen(cmd_asm, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    
    try:
        while process.poll() is None:
            current_mem = get_process_memory(process.pid)
            if current_mem > max_mem:
                max_mem = current_mem
            time.sleep(0.02)

    except Exception:
        process.kill()
        return 0.0, 0.0, 0.0

    duration = time.time() - start_time

    if process.returncode != 0:
        stderr_out = process.stderr.read()
        print(f"\n[ERROR] Assembly go wrong for {input_file}:\n{stderr_out}")
        return 0.0, duration, max_mem

    cmd_eval = ["bash", "evaluate.sh", output_file]
    result = subprocess.run(cmd_eval, capture_output=True, text=True, errors='ignore')
    
    score = 0.0
    match = re.search(r"Łączna ocena:\s+([\d\.]+)", result.stdout)
    if match:
        score = float(match.group(1))
        
    # Clear up
    if os.path.exists(output_file):
        os.remove(output_file)
        
    return score, duration, max_mem


def main():
    if not os.path.exists("test_manifest.csv"):
        return None

    print(f"{'Cov (x)':<7} | {'Err (%)':<7} | {'Score':<7} | {'Base':<5} | {'Diff':<7} | {'Time (s)':<8} | {'RAM (MB)':<8} | {'Status'}")
    print("=" * 85)
    
    results_diff = []
    wins = 0
    total = 0

    with open("test_manifest.csv", "r") as f:
        header = next(f, None)
        if not header:
            return None

        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split(',')
            if len(parts) < 3:
                continue
            
            filename, err_str, cov_str = parts
            
            err_float = float(err_str)
            cov_float = float(cov_str)
            err_percent = err_float * 100
            
            baseline = get_baseline(err_percent)
            
            try:
                score, duration, ram = run_assembly_and_eval(filename)
            except Exception as e:
                print(f"Exception processing {filename}: {e}")
                score, duration, ram = 0.0, 0.0, 0.0
            
            diff = score - baseline
            
            # Status
            if diff >= 0:
                status = "PASS"
                wins += 1
            else:
                status = "FAIL"
            
            diff_str = f"{diff:+.3f}"
            print(f"{cov_float:<7.1f} | {err_percent:<7.1f} | {score:<7.3f} | {baseline:<5.2f} | {diff_str:<7} | {duration:<8.2f} | {ram:<8.2f} | {status}")
            results_diff.append(diff)
            total += 1

    print("=" * 85)
    print(f"SUMMARY:")
    if total > 0:
        print(f"Beat baseline in {wins}/{total} trials ({wins/total*100:.1f}%).")
        print(f"Average score improvement: {sum(results_diff)/len(results_diff):+.4f}")


if __name__ == "__main__":
    main()
