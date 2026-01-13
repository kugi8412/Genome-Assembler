#!/usr/bin/env python3
# grid_search.py
# -*- coding: utf-8 -*-

import os
import re
import itertools
import subprocess

# Config grid search
GRID = {
    "k": [17, 19, 21, 23],
    "min_cov": [1, 2],
    "beam_width": [20, 50, 100]
}
OUTPUT_DIR = "grid_results"
READS_DIR = "training/reads"
DATASETS = ["reads1.fasta", "reads2.fasta", "reads3.fasta"]


def run_assembly_forced(reads_file, output_file, k, cov, beam):
    cmd = [
        "python3", "assembly.py",
        reads_file,
        output_file,
        "--force_k", str(k),
        "--force_min_cov", str(cov),
        "--beam_width", str(beam),
        "--min_contig_len", "300"
    ]

    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_evaluation(contigs_file):
    cmd = ["bash", "evaluate.sh", contigs_file]
    result = subprocess.run(cmd, capture_output=True, text=True, errors='ignore')
    match = re.search(r"Łączna ocena:\s+([\d\.]+)", result.stdout)
    if match:
        return float(match.group(1))
    return 0.0

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"{'Dataset':<15} | {'K':<3} | {'Cov':<3} | {'Beam':<4} | {'Score':<8}")
    print("-" * 55)

    for k, cov, beam in itertools.product(GRID["k"], GRID["min_cov"], GRID["beam_width"]):
        for dataset in DATASETS:
            input_path = os.path.join(READS_DIR, dataset)
            output_file = os.path.join(OUTPUT_DIR, f"temp_{dataset}")
            
            try:
                run_assembly_forced(input_path, output_file, k, cov, beam)
                score = run_evaluation(output_file)
                print(f"{dataset:<15} | {k:<3} | {cov:<3} | {beam:<4} | {score:.4f}")
                
            except subprocess.CalledProcessError:
                print(f"{dataset:<15} | {k:<3} | {cov:<3} | {beam:<4} | ERROR")
            except Exception as e:
                print(f"Error: {e}")

    print("\nThe grid search is complete. Check the results above.")

if __name__ == "__main__":
    main()
