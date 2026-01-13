#!/usr/bin/env python3
# generate_data.py
# -*- coding: utf-8 -*-


import os
import sys
import random
import numpy as np

# Config
REF_FILE = "training/reference/reference.fasta"
OUTPUT_DIR = "stress_test_data"
READ_LEN = 80


def load_reference(filepath: str) -> str:
    if not os.path.exists(filepath):
        print(f"[ERROR]: Not find file: {filepath}")
        sys.exit(1)

    with open(filepath, 'r') as f:
        seq = ""
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()

    return seq

def mutate_sequence(seq: str, error_rate: float) -> str:
    """Make SNP mutations in the sequence at the given error rate.
    """
    if error_rate <= 0:
        return seq
    
    bases = ['A', 'C', 'G', 'T']
    read_list = list(seq)
    num_errors = int(len(seq) * error_rate)

    if num_errors == 0 and random.random() < (len(seq) * error_rate):
        num_errors = 1

    if num_errors > 0:
        indices = random.sample(range(len(seq)), num_errors)
        for idx in indices:
            original_base = read_list[idx]
            options = [b for b in bases if b != original_base]
            read_list[idx] = random.choice(options)
            
    return "".join(read_list)

def generate_dataset(ref_seq: str,
                     target_cov: float,
                     target_err: float,
                     filename: str
                     ) -> int:
    ref_len = len(ref_seq)

    # Coverage = (NumReads * ReadLen) / GenomeLen
    # NumReads = (Coverage * GenomeLen) / ReadLen
    num_reads = int((target_cov * ref_len) / READ_LEN)
    reads_data = []
    
    for i in range(num_reads):
        start_pos = random.randint(0, ref_len - READ_LEN)
        fragment = ref_seq[start_pos : start_pos + READ_LEN]
        read_seq = mutate_sequence(fragment, target_err)
        reads_data.append(read_seq)
    
    # Shuffle reads
    random.shuffle(reads_data)
    
    with open(filename, 'w') as f:
        for i, seq in enumerate(reads_data):
            f.write(f">read_{i}\n{seq}\n")
            
    return num_reads

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    manifest = []  
    ref_seq = load_reference(REF_FILE)  
    coverage_steps = np.arange(5.0, 10.1, 0.5)

    # reads1 (circa 1%); reads2 (1-3%); reads3 (3-5%) + skrajne przypadki
    error_points = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.055, 0.06]
    count = 0
    for cov in coverage_steps:
        for err in error_points:
            filename = os.path.join(OUTPUT_DIR, f"sim_cov{cov:.1f}_err{err*100:.1f}_{count}.fasta")
            real_reads_count = generate_dataset(ref_seq, cov, err, filename)
            manifest.append((filename, err, cov))
            sys.stdout.write(f"\rGenerate: {count+1}/{len(coverage_steps)*len(error_points)} -> Cov={cov:.1f}x, Err={err*100:.1f}%, Reads={real_reads_count}")
            sys.stdout.flush()
            count += 1

    with open("test_manifest.csv", "w") as f:
        f.write("filename,error_rate,coverage\n")
        for fname, err, cov in manifest:
            f.write(f"{fname},{err},{cov}\n")
            
if __name__ == "__main__":
    main()
