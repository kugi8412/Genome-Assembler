#!/usr/bin/env python3
# assembly.py
# -*- coding: utf-8 -*-

import sys
import argparse
import statistics
import collections


# Configs
BEAM_WIDTH_DEFAULT = 20  
MIN_CONTIG_LEN_DEFAULT = 300 


class Assembler:
    def __init__(self,
                 beam_width: int,
                 min_contig_len: int
                 ) -> None:
        self.reads = []
        self.graph = collections.defaultdict(dict)
        self.rev_graph = collections.defaultdict(dict)
        self.out_degree = collections.defaultdict(int)
        self.in_degree = collections.defaultdict(int)

        self.k = 21 
        self.min_cov = 2
        self.min_correction_cov = 3
        self.beam_width = beam_width
        self.min_contig_len = min_contig_len
        self.skip_tip_removal = False

    def load_reads(self, filepath):
        with open(filepath, 'r') as f:
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq: self.reads.append("".join(seq))
                    seq = []
                else:
                    seq.append(line)

            if seq: self.reads.append("".join(seq))

    def analyze_histogram(self):
        """Analyses the k-mer histogram for K=15.
        """
        sample_k = 15
        counts = collections.defaultdict(int)
        for read in self.reads:
            if len(read) < sample_k:
                continue

            for i in range(len(read) - sample_k + 1):
                kmer = read[i:i+sample_k]
                counts[kmer] += 1
        
        if not counts:
            return 0, 0

        # Singletons
        singletons = sum(1 for v in counts.values() if v == 1)
        total_unique = len(counts)
        
        # Error estimation
        error_ratio = singletons / total_unique if total_unique > 0 else 0
        
        # Coverage estimation
        solid_counts = [v for v in counts.values() if v > 1]
        est_cov = statistics.median(solid_counts) if solid_counts else 0.0
        
        return est_cov, error_ratio

    def auto_tune_parameters(self):
        """Selecting parameters based on histogram analysis.
        """
        est_cov, error_ratio = self.analyze_histogram()
    
        if est_cov >= 7.5:
            self.k = 21
            self.min_cov = 2
            self.min_correction_cov = 3
            self.beam_width = 20
            self.skip_tip_removal = False
        elif 4.0 <= est_cov < 7.5:
            if error_ratio < 0.25:
                self.k = 21
                self.min_cov = 1
                self.min_correction_cov = 1
                self.beam_width = 50
                self.skip_tip_removal = False
            else:
                self.k = 17
                self.min_cov = 2
                self.min_correction_cov = 2
                self.beam_width = 60
                self.skip_tip_removal = False
        else:
             self.k = 15
             self.min_cov = 1


    def count_kmers(self, reads_source: list):
        counts = collections.defaultdict(int)
        for read in reads_source:
            if len(read) < self.k:
                continue

            for i in range(len(read) - self.k + 1):
                kmer = read[i:i+self.k]
                counts[kmer] += 1

        return counts

    def error_correction(self):
        if self.min_correction_cov <= 1:
            return None

        raw_counts = self.count_kmers(self.reads)
        correction_thresh = self.min_correction_cov
        solid_kmers = {k for k, v in raw_counts.items() if v >= correction_thresh}
        
        corrected_reads = []
        bases = ['A', 'C', 'G', 'T']
        
        for read in self.reads:
            read_list = list(read)
            for i in range(len(read) - self.k + 1):
                kmer = "".join(read_list[i:i+self.k])
                if kmer not in solid_kmers:
                    candidates = []
                    for j in range(self.k):
                        orig = kmer[j]
                        for b in bases:
                            if b == orig:
                                continue

                            cand = kmer[:j] + b + kmer[j + 1:]
                            if cand in solid_kmers:
                                candidates.append((j, b, raw_counts[cand]))

                    if candidates:
                        best_cand = max(candidates, key=lambda x: x[2])
                        read_list[i + best_cand[0]] = best_cand[1]

            corrected_reads.append("".join(read_list))

        self.reads = corrected_reads

    def build_graph(self):
        self.graph.clear()
        self.rev_graph.clear()
        self.out_degree.clear()
        self.in_degree.clear()

        final_counts = self.count_kmers(self.reads)

        for read in self.reads:
            for i in range(len(read) - self.k):
                k1 = read[i:i+self.k]
                k2 = read[i+1: i+1+self.k]
                
                if final_counts[k1] < self.min_cov or final_counts[k2] < self.min_cov:
                    continue
                
                if k2 not in self.graph[k1]:
                    self.graph[k1][k2] = 0

                self.graph[k1][k2] += 1
                
                if k1 not in self.rev_graph[k2]:
                    self.rev_graph[k2][k1] = 0

                self.rev_graph[k2][k1] += 1

        for u, neighbors in self.graph.items():
            self.out_degree[u] = len(neighbors)
        for v, parents in self.rev_graph.items():
            self.in_degree[v] = len(parents)

    def remove_edge(self, u, v):
        if u in self.graph and v in self.graph[u]:
            del self.graph[u][v]
            if not self.graph[u]:
                del self.graph[u]

            self.out_degree[u] -= 1
        
        if v in self.rev_graph and u in self.rev_graph[v]:
            del self.rev_graph[v][u]
            if not self.rev_graph[v]:
                del self.rev_graph[v]

            self.in_degree[v] -= 1

    def remove_tips(self):
        changed = True
        while changed:
            changed = False
            candidates = [n for n in self.graph if self.out_degree[n] == 0] + \
                         [n for n in self.rev_graph if self.in_degree[n] == 0]
            seen = set()
            for node in candidates:
                if node in seen:
                    continue

                seen.add(node)
                
                if self.out_degree[node] == 0 and self.in_degree[node] > 0:
                    parents = list(self.rev_graph.get(node, {}).items())
                    for parent, weight in parents:
                        is_weak = weight < self.min_cov 
                        has_better_alternative = False
                        if self.out_degree[parent] > 1:
                            for alt_child, alt_weight in self.graph[parent].items():
                                if alt_child != node and alt_weight >= weight:
                                    has_better_alternative = True
                                    break

                        if is_weak or has_better_alternative:
                            self.remove_edge(parent, node)
                            changed = True
 
                elif self.in_degree[node] == 0 and self.out_degree[node] > 0:
                    children = list(self.graph.get(node, {}).items())
                    for child, weight in children:
                        is_weak = weight < self.min_cov
                        has_better_alternative = False
                        if self.in_degree[child] > 1:
                            for alt_parent, alt_weight in self.rev_graph[child].items():
                                if alt_parent != node and alt_weight >= weight:
                                    has_better_alternative = True
                                    break
    
                        if is_weak or has_better_alternative:
                            self.remove_edge(node, child)
                            changed = True

    def simplify_bubbles(self):
        for u in list(self.graph.keys()):
            if self.out_degree[u] < 2:
                continue

            neighbors = self.graph[u]
            paths = [] 
            for v, w_uv in neighbors.items():
                if v in self.graph:
                    for w in self.graph[v]:
                         paths.append({'v': v, 'weight': w_uv, 'w': w})

            targets = collections.defaultdict(list)

            for p in paths:
                targets[p['w']].append(p)

            for w, group in targets.items():
                if len(group) > 1:
                    best_path = max(group, key=lambda x: x['weight'])
                    for p in group:
                        if p['v'] != best_path['v']:
                            self.remove_edge(u, p['v'])

    def get_contigs_beam_search(self):
        contigs = []
        global_visited_edges = set() 
        nodes = list(self.graph.keys())
        
        start_nodes = sorted(nodes, key=lambda n: (self.in_degree[n] == 0,
                                                   sum(self.graph[n].values())),
                                                   reverse=True
                    )

        processed_starts = set()

        for start_node in start_nodes:
            if start_node in processed_starts:
                continue

            if any((u, start_node) in global_visited_edges for u in self.rev_graph.get(start_node, {})):
                 continue

            initial_beam = (0, start_node, [start_node], set())
            active_beams = [initial_beam]
            best_finished_path = None 
            
            while active_beams:
                next_candidates = []
                made_progress = False
                for score, curr, path, local_visited in active_beams:
                    neighbors = self.graph.get(curr, {})
                    if not neighbors:
                        final_score = score + len(path) * 0.1
                        if best_finished_path is None or final_score > best_finished_path[0]:
                            best_finished_path = (final_score, path)
                        continue

                    for neighbor, weight in neighbors.items():
                        edge = (curr, neighbor)
                        if edge in global_visited_edges:
                            continue

                        if edge in local_visited:
                            continue
                        
                        made_progress = True
                        new_score = score + weight
                        if self.out_degree[neighbor] == 0:
                            new_score -= 5
                        else:
                            new_score += 1

                        if self.in_degree[neighbor] > 1:
                            new_score -= 0.5 
                            
                        new_visited = local_visited.copy()
                        new_visited.add(edge)
                        next_candidates.append((new_score, neighbor, path + [neighbor], new_visited))
                
                if not made_progress:
                    curr_best = max(active_beams, key=lambda x: x[0])
                    curr_score = curr_best[0] + len(curr_best[2]) * 0.1
                    if best_finished_path is None or curr_score > best_finished_path[0]:
                        best_finished_path = (curr_score, curr_best[2])
                    
                next_candidates.sort(key=lambda x: x[0], reverse=True)
                active_beams = next_candidates[:self.beam_width]

            if best_finished_path:
                final_score, final_path_nodes = best_finished_path
                contig_seq = final_path_nodes[0]
                for i in range(1, len(final_path_nodes)):
                    contig_seq += final_path_nodes[i][-1]
                
                if len(contig_seq) >= self.min_contig_len:
                    contigs.append(contig_seq)
                    for i in range(len(final_path_nodes) - 1):
                        global_visited_edges.add((final_path_nodes[i], final_path_nodes[i + 1]))
                    for node in final_path_nodes:
                        processed_starts.add(node)

        return contigs


def main():
    parser = argparse.ArgumentParser(description="Genome Assembler")
    parser.add_argument("input_file", type=str, help="Input FASTA file")
    parser.add_argument("output_file", type=str, help="Output FASTA file")
    parser.add_argument("--beam_width", type=int, default=BEAM_WIDTH_DEFAULT)
    parser.add_argument("--min_contig_len", type=int, default=MIN_CONTIG_LEN_DEFAULT)
    parser.add_argument("--force_k", type=int, default=None)
    parser.add_argument("--force_min_cov", type=int, default=None)
    parser.add_argument("--force_min_corr", type=int, default=None)

    args = parser.parse_args()

    assembler = Assembler(beam_width=args.beam_width, min_contig_len=args.min_contig_len)
    
    try:
        assembler.load_reads(args.input_file)
    except FileNotFoundError:
        print(f"Error: File {args.input_file} not found.")
        sys.exit(1)
    
    assembler.auto_tune_parameters()

    if args.force_k:
        assembler.k = args.force_k
    if args.force_min_cov:
        assembler.min_cov = args.force_min_cov
    if args.force_min_corr:
        assembler.min_correction_cov = args.force_min_corr

    assembler.error_correction()
    assembler.build_graph()
    
    if not assembler.skip_tip_removal:
        assembler.remove_tips()
        
    assembler.simplify_bubbles()
    contigs = assembler.get_contigs_beam_search()
    
    with open(args.output_file, 'w') as f:
        count = 1
        for c in contigs:
            if len(c) >= args.min_contig_len:
                f.write(f">contig_{count}\n{c}\n")
                count += 1

if __name__ == "__main__":
    main()

