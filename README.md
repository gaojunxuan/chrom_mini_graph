# chrom_mini_graph 

chrom_mini_graph creates a chromatic (coloured) minimizer pangenome graph. 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.

### Install

```
git clone https://github.com/bluenote-1577/chrom_mini_graph
cd chrom_mini_graph
cargo build --release
./target/release/chrom_mini_graph generate test_refs/*
./target/release/chrom_mini_graph map serialized_mini_graph.bin test_reads/hg_01243_pacbio_reads.fastq -b hg-01243-reads.bam > output.txt
```

1. `cargo build --release` first builds the **chrom_mini_graph** binary, which is found in the ./target/release/ directory.
2. The `chrom_mini_graph generate` command generates a coloured minimizer pangenome graph. 
3. The `chrom_mini_graph map` command chains onto the output graph and produces an alignment.
4. The resulting chain is used for alignment and is output to `hg-01243-reads.bam`. 

* 6 reference 1M bp segments of chromosome 20 are provided in the test_ref folder. 
* Simulated PacBio CLR reads for hg01243 are available in the test_reads folder. 

# Using chrom_mini_graph

## generate

`chrom_mini_graph generate ref_1.fasta ref_2.fasta ... -o output_from_generate` to create a coloured minimizer pangenome graph for references ref_1.fasta, ref_2.fasta, etc. 

* The window size can be easily modified in the `src/bin/chrom_mini_graph.rs` file. The default value is 16.
* Outputs a \*.bin file to be used for mapping and other auxillary information; see below. 

## map

A proof of concept read-to-graph chainer by chaining minimizers in the read onto the graph without knowledge of colour and then outputting the scores corresponding to each colouring.

`chrom_mini_graph map output_from_generate.bin your_reads.fastq -b bam_name.bam > output.txt`

* In stdout, information about each alignment is output.
* Chaining scores are output in the form (Colour Bits, Score , X, X) where the colour bits encode haplotypes. For example:  [(32, 3727.0, 0), (16, 3759.0, 0), (8, 3724.0, 0), (4, 3758.0, 0), (2, 3949.0, 0), (1, 3821.0, 0)]
indicates that the colour corrresponding to 2 has the highest score; this is expected since hg01243 corresponds to the second smallest bit. 
* For the best haplotype (if there is more than one, we pick one at random) according to chaining, we align to the haplotype using [block-aligner](https://github.com/Daniel-Liu-c0deb0t/block-aligner) and output results in the bam file.

## Issues

1. Aligning two big contigs (chromosomes) takes a long time right now; > 2000 seconds. This is caused by extremely repetitive kmers creating too many anchors. Will work on removing repetitive k-mers later on (i.e. masking repetitive kmers).
2. If two contigs are dissimilar, the graph generation may be very poor. We find the best alignment and align no matter what; we don't check if the alignment is actually good or not. 
3. I have not optimized the alignment. It seems like it is more likely to fail in noisy regions than minimap2 alignment. 
4. Only one alignment is output per read. No supplementary/secondary alignments are output.

## Outputs from `generate`

### full_mini_graph.csv 

This is a CSV file which indicates the edges for the entire graph. 

### simpilified_mini_graph.csv

This is a simplified version of the full_mini_graph.csv file. Linear paths are shortened. 

### serialized_mini_graph.json

This is the JSON serialization for the entire graph. See [this document](https://docs.google.com/document/d/1oRHjPgP-Bh9UkySCduWIl5yCpfiLVEoSnRdzdx4a7-Y/edit?usp=sharing) for how to deserialize the graph. **IMPORTANT:** For colouring, the most significant bit corresponds to the first genome in the command, and the least significant bit corresponds to the last genome in the command. 

### serialized_mini_graph.bin

Binary serializeation of the graph used in the `map` subcommand.

## Outputs from `map`

### read_anchor_hits.txt

Gives the ids of the nodes corresponding to the best chain for each read. 

### Visualization

To quickly visualize the graph, do `python visualize_graph.py full_mini_graph.csv 26000 27000`. All nodes with ids between 26000 and 27000 will be visualized here. You need networkx, graphviz, matplotlib installed.
