import json
import sys, getopt

def main(argv):
    opts, args = getopt.getopt(argv, "khe:i:o:",["input=","output=", "edgelst="])
    in_file = ""
    out_file = ""
    edgelst_file = ""
    k = 16
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: convert_to_gfa -i input.json -o output.gfa [-e output_edgelist.txt -k kmer_len]")
            return
        elif opt == '-i' or opt == '--input':
            in_file = arg
        elif opt == '-o' or opt == '--output':
            out_file = arg
        elif opt == '-e' or opt == '--edgelst':
            edgelst_file = arg
        elif opt == '-k':
            if arg.isdigit():
                k = int(arg)
            else:
                print("[Error] kmer length has to be an integer")
                return
    if in_file == "" or out_file == "":
        print("Usage: convert_to_gfa -i input.json -o output.gfa [-e output_edgelist.txt]")
        return
        
    with open(in_file) as f:
        data = json.load(f)

    if data:
        dim = len(data)
        print("[Info] dim =", dim)
        edge_list = []
        node_lengths = []
        for i, node in enumerate(data):
            node_lengths.append(len(node['internal_nodes']) * k)
            if node['child_nodes']:
                for j, child_idx in enumerate(node['child_nodes']):
                    dist = node['child_edge_distance'][j][0]
                    edge_list.append((node['pseudo_id'], child_idx, dist))
        if edgelst_file != "": 
            with open(edgelst_file, 'a') as ef:
                for u, v, dist in edge_list:
                    ef.write('{},{},{}\n'.format(u, v, dist))
        num_nodes = len(data)
        num_edges = len(edge_list)
        print("[Info] n = {}, m = {}".format(num_nodes, num_edges))
        with open(out_file, 'w') as gfa_file:
            gfa_file.write('H\tVN:Z:2.0\n')
            for i in range(num_nodes):
                gfa_file.write('S\t{}\t*\tLN:i:{}\n'.format(i, node_lengths[i]))
            for j in range(num_edges):
                edge = edge_list[j]
                gfa_file.write('L\t{}\t+\t{}\t+\t0M\n'.format(edge[0], edge[1], 1, 1))
    else:
        print("[Error] Failed to load serialized json file.")
        return

if __name__ == "__main__":
    main(sys.argv[1:])
