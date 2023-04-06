#!/usr/bin/python3

import subprocess

for i in range(2,201,2):
    name = f'data/exact_{i:03}.gr'
    with open(name,'r') as f:
        _,_,n,m = f.readline().split()
    print(i//2, name, n,m)
    try:
        subprocess.run(['./cmake-build-release/source/example'], stdin=open(name), timeout=60)
    except Exception:
        print('timeout')
    
