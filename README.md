# chrom_mini_graph 

chrom_mini_graph creates a chromatic (coloured) minimizer pangenome graph. 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.

### Install

```
git clone https://github.com/bluenote-1577/chrom_mini_graph
cd chrom_mini_graph
cargo build --release
./target/release/chrom_mini_graph test_ref/*
```

`cargo build --release` builds the **chrom_mini_graph** binary, which is found in the ./target/release/ directory. 

## Using chrom_mini_graph

`./target/release/chrom_mini_graph ref_1.fasta ref_2.fasta ...` to create a coloured minimizer pangenome graph for references ref_1.fasta, ref_2.fasta, etc. 

6 reference 1M bp segments of chromosome 20 are provided in the test_ref folder. 

To modify the parameters such as w, k, or even use syncmers instead of minimizers, see the parameters in the `src/bin/chrom_mini_graph.rs` file.

## Output

### full_mini_graph.csv 

This is a CSV file which indicates the edges for the entire graph. 

### simpilified_mini_graph.csv

TODO: This is a simplified version of the full_mini_graph.csv file. Not working for multi-edges right now, I will fix this later. 

### serialized_mini_graph.json

This is the JSON serialization for the entire graph. See [this document](https://docs.google.com/document/d/1oRHjPgP-Bh9UkySCduWIl5yCpfiLVEoSnRdzdx4a7-Y/edit?usp=sharing) for how to deserialize the graph.

### Visualization

To quickly visualize the graph, do `python visualize_graph.py full_mini_graph.csv`. The parameters in the `visualize_graph.py` file can be changed to visualize different sections of the graph.
