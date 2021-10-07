import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph
import sys
import graphviz

use_graphviz = False
a1 = 26000
a2 = 27000
b1 = 26000
b2 = 27000

if not use_graphviz:
    file = open(sys.argv[1],'r')
    G = nx.Graph()
    for line in file:
        sp = line.split(',')
        n1 = int(sp[0].split('-')[0])
        n2 = int(sp[1].split('-')[0])
        #w = int(sp[2])
        if a1 < n1 < a2 and b1< n2 < b2:
            #print(n1,n2,w)
            #G.add_edge(n1,n2,weight=w)
            G.add_edge(n1,n2)

    print(G.edges)
    nx.draw(G, node_size=5, pos = nx.nx_pydot.graphviz_layout(G))
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
        if a1 < int(n1) < a2 and b1 < int(n2) < b2:
            G.edge(n1,n2, label=w)
    G.render()






