#!/bin/bash
#$ -P simpsonlab
#$ -N chrom_minigraph_benchmark_generate
#$ -cwd
#$ -l h_vmem=32G
#$ -m eas
#$ -M your_email_address

# Load shell
source /oicr/local/Modules/default/init/sh

# Define variables
CMG_EXE=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
CHR11_FASTA=/u/jgao/jgao/projects/chrom_minigraph/chr11/*
CHR20_FASTA=/u/jgao/jgao/projects/chrom_minigraph/chr20/* 

# Run
$CMG_EXE generate $CHR11_FASTA