#!/usr/bin/python3

import subprocess

for i in range(2,101,2):
    name = f'../data/exact_{i:03}.gr'
    print(name)
    with open(name,'r') as f:
        _,_,n,m = f.readline().split()
        edges = list(map(lambda w: tuple(map(int,w.split())), f))
    s = '(' +  str(list(range(int(n)))) +', [{' + '}, {'.join(f'{u-1},{v-1}' for u,v in edges) + '}])'
    
    with open('tmp.dot', 'w') as f:
        subprocess.run(['../../modular-decomposition/build/mod_dec', s], stdout=f)

    trivial = True
    with open('tmp.dot', 'r') as f:
        md = f.read()
        if 'Series' in md or 'Parallel' in md:
            trivial = False
    
    if not trivial:
        subprocess.run(['dot', '-Tpng', 'tmp.dot', '-o', f'{i}.png'])
            
