#$ -N chrom_minigraph_benchmark
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=144000
#$ -m eas
#$ -M your_email_address

# For each ref length
# 	create reference graph using 10 references
# 	report time and mem usage of graph generation time
# 	for each of the 10 composing read of the reference
# 		map read to the graph
# 		report time and mem usage for mapping
# 		sort bam file
# 		report coverage and depth

CMG_FOLDER=/u/jgao/jgao/code/chrom_mini_graph
CMG=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
GENERATED_FILE=/u/jgao/jgao/code/chrom_mini_graph/serialized_mini_graph.bin
BAM_FILE=/u/jgao/jgao/code/chrom_mini_graph/output.bam
SORTED_BAM=/u/jgao/jgao/code/chrom_mini_graph/sorted_output.bam
SIZES=(100000, 250000, 500000, 750000, 1000000, 2500000, 5000000, 7500000, 10000000)
PYTHON3=/u/jgao/miniconda3/bin/python3

module load samtools/1.17

for i in ${SIZES[@]}
do
    qsub -N "chrom_minigraph_generate_size_${i}" -cwd -l h_vmem=32G -l h_rt=72000 -m eas -M your_email_address -q u20.q $CMG generate /u/jgao/jgao/code/rand_seq/${i}/refs/*.fasta
    for j in {0..9}
    do
        qsub -N "chrom_minigraph_map_size_${i}_ref_${j}" -cwd -l h_vmem=32G -l h_rt=72000 -m eas -M your_email_address -q u20.q $CMG map -a --trace $GENERATED_FILE /u/jgao/jgao/code/rand_seq/${i}/reads/ref${j}_len${i}.fasta.fastq
        samtools sort $BAM_FILE -o $SORTED_BAM
        samtools coverage $SORTED_BAM > $CMG_FOLDER/coverage_report_size_${i}_ref_${j}.txt
        $PYTHON3 $CMG_FOLDER/analyze_mapping_accuracy.py ref${j}_len${i} $CMG_FOLDER/best_genome_reads.txt
    done
done

