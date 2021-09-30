import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph
import sys
import graphviz

use_graphviz = True

if not use_graphviz:
    file = open(sys.argv[1],'r')
    G = nx.MultiDiGraph()
    for line in file:
        sp = line.split(',')
        n1 = int(sp[0])
        n2 = int(sp[1])
        w = int(sp[2])
        if n1 < 10000 and n2 < 10000:
            print(n1,n2,w)
            G.add_edge(n1,n2,weight=w)

    print(G.edges)
    #nx.draw(G, node_size=5, pos = nx.nx_pydot.graphviz_layout(G),connectionstyle='arc3, rad = 0.1')
    #plt.show()
    A = to_agraph(G)
    A.layout('dot')
    A.draw('multi.png')
else:
    file = open(sys.argv[1],'r')
    G = graphviz.Digraph(format='png', strict=False)
    for line in file:
        sp = line.split(',')
        n1 = sp[0]
        n2 = sp[1]
        w = sp[2]
        print(n1,n2,w)
        if int(n1) < 10000 and int(n2) < 10000:
            G.edge(n1,n2, label=w)

    G.render()






