#!/usr/bin/python3

import subprocess
import sys
import shutil

for num in range(1,101):
    name = f'data/exact_{2*num:03}.gr'

    print(num, name, flush=True)
    res = subprocess.run(['./cmake-build-release/main', '-v'], 
            stdin=open(name), 
            stdout=subprocess.DEVNULL)

    shutil.move('graph.dot', f'tmp/{num}.dot')

