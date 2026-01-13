#!/usr/bin/env python3
# benchmark.py
# -*- coding: utf-8 -*-


import os
import sys
import re
import time
import subprocess


PARAMS = {
    "--beam_width": "20",
    "--min_contig_len": "300"
}
OUTPUT_DIR = "final_outs"
READS_DIR = "training/reads"
DATASETS = ["reads1.fasta", "reads2.fasta", "reads3.fasta"]
# Resource example limits
LIMIT_RAM_MB = 500.0
LIMIT_TIME_SEC = 3600.0


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


def run_evaluation(contigs_file):
    cmd = ["bash", "evaluate.sh", contigs_file]
    result = subprocess.run(cmd, capture_output=True, text=True, errors='ignore')
    match = re.search(r"Łączna ocena:\s+([\d\.]+)", result.stdout)
    if match:
        return float(match.group(1))

    return 0.0


def run_monitored_assembly(input_file, output_file, custom_params=None):
    """ Run assembly.py with resource monitoring.
    """
    cmd = ["python3", "assembly.py", input_file, output_file]
    params_to_use = custom_params if custom_params else PARAMS
    for key, value in params_to_use.items():
        cmd.extend([key, str(value)])

    start_time = time.time()
    process = subprocess.Popen(cmd)
    max_mem = 0.0
    
    try:
        while process.poll() is None:
            current_mem = get_process_memory(process.pid)
            if current_mem > max_mem:
                max_mem = current_mem
            time.sleep(0.02)
            
            if (time.time() - start_time) > LIMIT_TIME_SEC + 60:
                process.kill()
                return max_mem, time.time() - start_time, False

    except KeyboardInterrupt:
        process.kill()
        sys.exit(1)

    duration = time.time() - start_time
    success = (process.returncode == 0)
    return max_mem, duration, success


def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    results = []

    print(f"{'='*60}")
    print(f"BENCHMARK TESTY ASSEMBLY.PY")
    print(f"{'='*60}\n")

    for dataset in DATASETS:
        input_path = os.path.join(READS_DIR, dataset)
        output_path = os.path.join(OUTPUT_DIR, f"contigs_{dataset}")
        
        if not os.path.exists(input_path):
            continue

        # 1. Uruchomienie Assemblera
        peak_mem, duration, success = run_monitored_assembly(input_path, output_path)
        
        # 2. Ocena (Evaluate)
        score = 0.0
        if success:
            score = run_evaluation(output_path)
        
        status_icon = "✅" if success else "❌"
        mem_ok = "OK" if peak_mem <= LIMIT_RAM_MB else "NOT OK"
        time_ok = "OK" if duration <= LIMIT_TIME_SEC else "NOT OK"
        
        print(f"   Results Assembly: {status_icon}")
        print(f"   📊  SCORE:       {score:.4f}")
        print(f"   ⏱️   TIME:        {duration:.2f} s  ({time_ok})")
        print(f"   💾 MEMORY:      {peak_mem:.2f} MB ({mem_ok})")
        print(f"{'-'*60}")
        results.append((dataset, peak_mem, duration, score))

    # Summary
    print("\n" + "="*70)
    print(f"{'Dataset':<15} | {'RAM (MB)':<10} | {'Time (s)':<10} | {'Score':<10}")
    print("-" * 70)
    
    total_time = 0
    max_peak_ram = 0
    total_score = 0
    
    for ds, mem, t, scr in results:
        print(f"{ds:<15} | {mem:<10.2f} | {t:<10.2f} | {scr:<10.4f}")
        total_time += t
        total_score += scr
        if mem > max_peak_ram: max_peak_ram = mem
        
    print("="*70)
    print(f"Sum of Time: {total_time:.2f} s")
    print(f"Max RAM:    {max_peak_ram:.2f} MB")
    print(f"Sum of Score:   {total_score:.4f}")

    if max_peak_ram < LIMIT_RAM_MB and total_time < LIMIT_TIME_SEC:
        print("\n✅ Benchmark succeeded.")
    else:
        print("\n⚠️  Benchmark failed resource limits.")


if __name__ == "__main__":
    main()
