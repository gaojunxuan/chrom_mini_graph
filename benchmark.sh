#!/bin/bash
#$ -P simpsonlab
#$ -N chrom_minigraph_benchmark
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=144000
#$ -m eas
#$ -M jgao@oicr.on.ca

CMG_FOLDER=/u/jgao/jgao/code/chrom_mini_graph
CMG=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
GENERATED_FILE=/u/jgao/jgao/code/chrom_mini_graph/serialized_mini_graph.bin
BAM_FILE=/u/jgao/jgao/code/chrom_mini_graph/output.bam
SORTED_BAM=/u/jgao/jgao/code/chrom_mini_graph/sorted_output.bam
#SIZES=(100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000)
SIZES=(100000 1000000)
PYTHON3=/u/jgao/miniconda3/bin/python3
SAMTOOLS=/.mounts/labs/simpsonlab/sw/samtools/1.17/bin/samtools
QSUB=/opt/uge-8.6/bin/lx-amd64/qsub
SIM_FOLDER=/u/jgao/jgao/projects/chrom_minigraph/sim

for i in ${SIZES[@]}
do
    $CMG_FOLDER/run_generate.sh "$SIM_FOLDER/${i}/refs/*.fasta" &> $CMG_FOLDER/log_generate_size_${i}.txt
    for j in {0..9}
    do
        $CMG_FOLDER/run_map.sh $GENERATED_FILE "$SIM_FOLDER/${i}/reads/ref${j}_len${i}.fasta.fastq" &> $CMG_FOLDER/log_map_size_${i}_ref_${j}.txt
        $SAMTOOLS sort $BAM_FILE -o $SORTED_BAM
        $SAMTOOLS coverage $SORTED_BAM > $CMG_FOLDER/coverage_report_size_${i}_ref_${j}.txt
        $PYTHON3 $CMG_FOLDER/analyze_mapping_accuracy.py ref${j}_len${i} $CMG_FOLDER/best_genome_reads.txt > accuracy_report_size_${i}_ref_${j}.txt
    done
done