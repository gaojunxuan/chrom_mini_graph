import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph
import sys
import graphviz
import json

use_graphviz = False

if not use_graphviz:
    file = open(sys.argv[1],'r')
    a1 = int(sys.argv[2])
    a2 = int(sys.argv[3])
    bit = 1
    json_string = ""
    graph_deserial = []
    reads = dict()
    color = False

    if len(sys.argv) > 4:
        color = True
        read_to_check = sys.argv[7]
        json_file = open(sys.argv[4],'r')
        read_anchor_file = open(sys.argv[6],'r')
        json_string = json_file.readline()
        graph_deserial = json.loads(json_string)
        bit = int(sys.argv[5])
        for read in read_anchor_file:
            splitted = read.split(':')
            spl = splitted[1].split(',')
            int_spl = [int(x) for x in spl];
            for node in spl:
                reads[splitted[0]] = int_spl


    G = nx.DiGraph()
    for line in file:
        sp = line.split(',')
        n1 = int(sp[0].split('-')[0])
        n2 = int(sp[1].split('-')[0])
        #w = int(sp[2])
        if a1 < n1 < a2 and a1< n2 < a2:
            #print(n1,n2,w)
            #G.add_edge(n1,n2,weight=w)
            G.add_edge(n1,n2)

    if color:
        order_to_id = dict()
        for node in graph_deserial:
            order_to_id[node['order']] = int(node['id'])
        cmap = []
        node_s = []
        for n_order in G:
            inside = False
            if order_to_id[n_order] in reads[read_to_check]:
                inside = True
                #print(n_order)
                node_s.append(105)
            else:
                node_s.append(5)

            node_color = graph_deserial[order_to_id[n_order]]['color']
            if node_color & bit == bit:
                cmap.append('red')
            else:
                cmap.append('blue')

        nx.draw(G, node_size=node_s, pos = nx.nx_pydot.graphviz_layout(G), node_color = cmap)
    else:
        #nx.draw(G, node_size=7, pos = nx.nx_pydot.graphviz_layout(G), with_labels=True)
        nx.draw(G, node_size=7, pos = nx.nx_pydot.graphviz_layout(G))
    plt.show()
else:
    file = open(sys.argv[1],'r')
    G = graphviz.Digraph(format='png', strict=False)
    for line in file:
        sp = line.split(',')
        n1 = sp[0]
        n2 = sp[1]
        w = sp[2]
        print(n1,n2,w)
        G.edge(n1,n2, label=w)
    G.render()






