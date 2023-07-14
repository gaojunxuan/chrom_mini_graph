import getopt
import sys


def main(argv):
    opts, args = getopt.getopt(argv, "ht", ["help", "count_ties"])
    count_ties = False
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Usage: analyze_mapping_accuracy.py [-t] <expected_chromosome> <best_genome_reads.txt>")
            return
        elif opt in ("-t", "--count_ties"):
            count_ties = True
        if len(args) != 2:
            print("Usage: analyze_mapping_accuracy.py [-t] <expected_chromosome> <best_genome_reads.txt>")
            return
    expected_chromosome = args[0]
    best_genome_reads_file = args[1]
    num_correct = 0
    num_total = 0
    print("Count ties: {}".format(count_ties))
    with open(best_genome_reads_file) as f:
        lines = f.read().split(">")[1:]
    for line in lines:
        line = line.split("\n")
        read_name = line[0]
        chrom_name_score_pairs = line[1:]
        mapping_pairs = []
        score_of_expected_best_mapping = 0
        for chrom_name_score_pair in chrom_name_score_pairs:
            if chrom_name_score_pair == "" or chrom_name_score_pair.isspace():
                continue
            chrom_name_score_pair = chrom_name_score_pair.split("\t")
            chrom_name = chrom_name_score_pair[0]
            score = chrom_name_score_pair[1]
            if chrom_name == expected_chromosome:
                score_of_expected_best_mapping = score
            mapping_pairs.append((chrom_name, score))
        # sort by score
        mapping_pairs.sort(key=lambda x: x[1], reverse=True)
        # get the best mapping
        best_mapping = mapping_pairs[0]
        # if top mapping is expected or tie for top mapping is expected
        if best_mapping[0] == expected_chromosome or \
                (count_ties and
                 best_mapping[1] == score_of_expected_best_mapping):
            num_correct += 1
        num_total += 1
    print("Accuracy: {}/{} = {}".format(num_correct,
                                        num_total,
                                        num_correct / num_total))


if __name__ == "__main__":
    main(sys.argv[1:])
