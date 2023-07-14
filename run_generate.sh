#!/bin/bash
#$ -P simpsonlab
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=144000
#$ -m eas
#$ -M jgao@oicr.on.ca

CMG_FOLDER=/u/jgao/jgao/code/chrom_mini_graph
CMG=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
GENERATED_FILE=/u/jgao/jgao/code/chrom_mini_graph/report/$1/serialized_mini_graph.bin
SIM_FOLDER=/u/jgao/jgao/projects/chrom_minigraph/sim

$CMG generate $SIM_FOLDER/$1/refs/*.fasta