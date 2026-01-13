# Genome Assemmbler
This repository presents implementations of an advanced short-read assembler based on de Bruijn graphs, equipped with an automatic parameter tuning mechanism.
The programme works by dynamically adjusting its strategy depending on the quality of the input data (coverage and error rate).

## Modules

### Auto-Tuning
Before starting the actual _de novo_ assembly, the program performs a quick scan of the reads using small k-mers of size 15 to build a frequency histogram.
Based on this, two key metrics are calculated:

- **EstCov**: The median number of k-mers occurring more than once.
- **ErrorRatio**: The ratio of unique k-mers to all unique sequences.

| Condition                         | Mode               | Parameters             | Goal                                                                  |
|----------------------------------|--------------------|------------------------|-----------------------------------------------------------------------|
| EstCov ≥ 7.5                       | STABLE HIGH        | K=21, MinCov=2         | High precision for dense data. Noise filtering.               |
| 4.0 ≤ EstCov < 7.5 (ErrorRatio < 25%) | PRECISION RESCUE   | K=21, MinCov=1         | Preserves sparse reads, large K for precision. |
| 4.0 ≤ Cov < 7.5 (ErrorRatio ≥ 25%)    | AGGRESSIVE FILTER  | K=17, MinCov=2         | Smaller K increases overlap probability, MinCov=2 removes noise.       |
| Cov < 4.0                      | DESPERATE          | K=15, MinCov=1         | Extremely low coverage. Desperate attempt to assemble anything.       |

### Error Correction

Performs spectral correction before graph construction (unless **Rescue** mode has disabled it):
- Identifies rare k-mers in reads (below the **MinCov** threshold).
- Searches for neighbours (Hamming distance = 1) that are stronger part of the genome.
- Replaces incorrect bases with correct ones, repairing reads _in silico_.

### Building de Bruijn graph
Constructs a directed multigraph where the nodes are [k-1]-mers and the edges are k-mers.
Filters edges with a weight below **MinCov** (unless we are in **Rescue** mode).
Creates a forward and backward graph structure for efficient navigation.
This approach consumes twice as much memory, but thanks to this,
every operation requiring backtracking or checking the parent is performed in constant time.

### Graph Cleaning
The raw graph contains artefacts resulting from sequencing errors. The programme uses two topological techniques:
- Remove Tips: Removes short branches ending in the nowhere, which are usually errors at the ends of reads.
- Simplifying Bubbles: Detects places where the path splits and immediately rejoins (polymorphisms or errors). It retains the path with higher coverage, removing the weaker branch.

### Contig Generation
Instead of a simple greedy approach, the assembler uses Beam Search with a width defined by beam width.
Starts with the nodes with the highest input/output degree (the most reliable points). Then, at each step, instead of analysing the single best path,
it analyses as many paths as the beam width. We use the following evaluation function for the assessment:

        +1 for extending by an edge.

        + Edge Weight (prefers paths with high coverage).

        -5 penalty for entering a Dead-End.

        -0.5 penalty for ambiguity (branching).

The heuristics used allow the algorithm to avoid local minima, which the classic greedy algorithm falls into.

### Usage
`python assembly.py [input_reads.fasta] [output.fasta]`

#### List of Parameters:

        input_file: Path to the FASTA file with reads.

        output_file: Path to save the resulting contigs.

        --beam_width: Beam width in Beam Search (default: 20).

        --min_contig_len: Minimum contig length in the result (default: 300).

        --force_k: Forces k-mer length (ignores auto-tuning).

        --force_min_cov: Forces the k-mer filtering threshold.

        --force_min_corr: Forces the error correction threshold.

#### Requirments
Scripts require `bowtie2` program and Python modules:
- pysam
- numpy

### Evaluation
`bash evaluate.sh [output.fasta]`
