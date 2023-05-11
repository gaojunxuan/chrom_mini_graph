import random
import sys, getopt

nucleotide = ['A', 'C', 'G', 'T']

def main(argv):
    opts, args = getopt.getopt(argv, "hl:i:d:v:f:o:",["length=","insert_length=", "deletion_length=", "inversion_length", "flanking_length", "seq_name="])
    out_file = ""
    length = 0
    insert_len = -1
    del_len = -1
    inv_len = -1
    flanking_len = -1
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: rand_seq -l length -i insert_len -d deletion_length -v inversion_length -f flanking_length -o output.txt")
            return
        elif opt == '-l':
            if arg.isdigit():
                length = int(arg)
            else:
                print("[Error] length has to be an integer")
                return
        elif opt == '-i':
            if arg.isdigit():
                insert_len = int(arg)
            else:
                print("[Error] insert_len has to be an integer")
                return
        elif opt == '-o':
            out_file = arg
        elif opt == '-d':
            if arg.isdigit():
                del_len = int(arg)
            else:
                print("[Error] deletion_length has to be an integer")
                return
        elif opt == '-v':
            if arg.isdigit():
                inv_len = int(arg)
            else:
                print("[Error] inversion_length has to be an integer")
                return
        elif opt == '-f':
            if arg.isdigit():
                flanking_len = int(arg)
            else:
                print("[Error] flanking_length has to be an integer")
                return
    if out_file == "":
        print("Usage: rand_seq -l length -i insert_len -d deletion_length -v inversion_length -f flanking_length -o output.txt")
        return
    if length == 0:
        print("[Error] length cannot be 0")
        return
    if insert_len == -1:
        print("[Error] specify insert_len")
        return
    if del_len == -1:
        print("[Error] specify deletion_length")
        return
    if inv_len == -1:
        print("[Error] specify inversion_length")
        return
    if flanking_len == -1:
        print("[Error] specify flanking_length")
        return
    
    
    seq = ''.join([random.choice(nucleotide) for _ in range(length)])
    insert = ''.join([random.choice(nucleotide) for _ in range(insert_len)])
    insert_pos = random.randint(0, length - insert_len)
    
    with open(out_file + ".fa", 'w') as f:
        # write base sequence
        f.write(">{}\n{}\n".format(out_file, seq))

    with open(out_file + "_with_insert.fa", 'w') as f:
        # write sequence with insert
        seq_insert = seq[:insert_pos] + insert + seq[insert_pos:]
        f.write(">{}_with_insert\n{}\n".format(out_file, seq_insert))

    # with open(out_file + "_with_deletion.fa", 'w') as f:
    #     del_len = random.randint(1, insert_len)
    #     del_pos = random.randint(0, length - del_len)
    #     seq_del = seq[:del_pos] + seq[del_pos + del_len:]
    #     f.write(">{}_with_deletion\n{}\n".format(out_file, seq_del))

    # with open(out_file + "_with_inversion.fa", 'w') as f:
    #     inv_len = random.randint(1, insert_len)
    #     inv_pos = random.randint(0, length - inv_len)
    #     seq_inv = seq[:inv_pos] + seq[inv_pos:inv_pos + inv_len][::-1] + seq[inv_pos + inv_len:]
    #     f.write(">{}_with_inversion\n{}\n".format(out_file, seq_inv))
    
    with open(out_file + ".fastq", 'w') as f:
        # write portion of base sequence
        start_pos = random.randint(0, length - insert_len)
        f.write("@{}\n{}\n+\n{}\n".format(out_file, seq[start_pos:start_pos + insert_len], '+' * insert_len))

        # write portion of the inserted sequence
        start_pos = random.randint(0, insert_len - 1)
        f.write("@{}_inserted_portion\n{}\n+\n{}\n".format(out_file, insert[start_pos:], '+' * insert_len))

        # write portion of the inserted sequence plus a flanking region
        insert_seq = insert[start_pos:]
        seq_insert = seq[:insert_pos] + insert + seq[insert_pos:]
        flank_seq = seq_insert[insert_pos + insert_len: insert_pos + insert_len + flanking_len]
        f.write("@{}_inserted_portion_with_flanking\n{}\n+\n{}\n".format(out_file, insert_seq + flank_seq, '+' * (insert_len + flanking_len)))

if __name__ == "__main__":
    main(sys.argv[1:])