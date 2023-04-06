#!/usr/bin/python3

import networkx as nx

for i in range(2,201,2):
    name = f'data/exact_{i:03}.gr'
    with open(name,'r') as f:
        _,_,n,m = f.readline().split()
        n = int(n)
        m = int(m)
        edges = list(map(lambda w: tuple(map(int,w.split())), f))
    print(name, n,m, round(2*m/n,2))

    G = nx.Graph(edges)
    CCs = list(nx.connected_components(G))
    print('comps: ', len(CCs), ' (max=', max(map(len,CCs)), ')')
    print()
    
