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
./target/release/chrom_mini_graph map serialized_mini_graph.bin test_reads/hg_01243_pacbio_reads.fastq > output.txt
```

`cargo build --release` builds the **chrom_mini_graph** binary, which is found in the ./target/release/ directory. 

## Using chrom_mini_graph

### generate

`./target/release/chrom_mini_graph generate ref_1.fasta ref_2.fasta ...` to create a coloured minimizer pangenome graph for references ref_1.fasta, ref_2.fasta, etc. 

6 reference 1M bp segments of chromosome 20 are provided in the test_ref folder. 

To modify the parameters such as w, k, or even use syncmers instead of minimizers, see the parameters in the `src/bin/chrom_mini_graph.rs` file.

### map

A proof of concept read-to-graph mapper by mapping reads onto the graph without knowledge of colour, and then outputting the scores corresponding to each colouring.

`./target/release/chrom_mini_graph map (output_from_generate.bin) (your_reads.fastq) > output.txt`

For each read, the scores corresponding to each color is output in the form (Score, Colour) where colour is an integer corresponding to a bit. The file read_anchor_hits.txt gives the chains for each read in terms of nodes (the numbers are _orders_ of the nodes, not the IDs) on the graph. 

## Issues

1. Aligning two big contigs (chromosomes) takes a long time right now; > 2000 seconds. This is caused by extremely repetitive kmers creating too many anchors. Will work on removing repetitive k-mers later on (i.e. masking repetitive kmers).
2. **Make sure there are no stretches of N's in the reference file!**. These Ns will cause chaining to take a long time, because the Ns are converted to As in the code and repetitive k-mers cause the number of anchors to blow up. I just removed the Ns from the reference.
3. If two contigs are dissimilar, the graph generation may be very poor. We find the best alignment and align no matter what; we don't check if the alignment is actually good or not. 
4. Reverse complements don't work in any shape or form so be careful.

## Outputs

### full_mini_graph.csv 

This is a CSV file which indicates the edges for the entire graph. 

### simpilified_mini_graph.csv

This is a simplified version of the full_mini_graph.csv file. Linear paths are shortened. 

### serialized_mini_graph.json

This is the JSON serialization for the entire graph. See [this document](https://docs.google.com/document/d/1oRHjPgP-Bh9UkySCduWIl5yCpfiLVEoSnRdzdx4a7-Y/edit?usp=sharing) for how to deserialize the graph. **IMPORTANT:** For colouring, the most significant bit corresponds to the first genome in the command, and the least significant bit corresponds to the last genome in the command. 

### Visualization

To quickly visualize the graph, do `python visualize_graph.py full_mini_graph.csv 26000 27000`. All nodes with ids between 26000 and 27000 will be visualized here. You need networkx, graphviz, matplotlib installed.
