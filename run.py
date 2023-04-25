#!/usr/bin/python3

import subprocess
import pandas as pd
import sys
import json
from datetime import datetime

date = datetime.today().strftime('%Y-%m-%d')
logfile = f'tmp/{date}.log'

df = pd.read_csv('results/current.csv', index_col='num')
newly_solved = []
newly_lost = []

for num in list(range(1,35))+[100]:
    name = f'data/exact_{2*num:03}.gr'

    prev_low = df.loc[num,'low']
    prev_sol = df.loc[num,'solved']
    prev_opt = df.loc[num,'opt']
    prev_space = df.loc[num,'space']
    prev_time = df.loc[num,'time']

    print(num, name, flush=True)
    timeout = False
    with open(logfile,'a') as f:
        f.write(f'\n{num}\n')
    try:
        res = subprocess.run(['./cmake-build-release/main', '-v'], 
            stdin=open(name), 
            stdout=subprocess.DEVNULL,
            stderr=open(logfile,'a'), 
            timeout=60)
    except Exception as e:
        print(e)
        timeout = True

    if not timeout:
        with open('run_info.json') as f:
            info = json.load(f)
        opt = info['width']
        spa = info['search_space']
        tme = info['running_time']
        low = info['lower_bound']
        if not prev_sol:
            print(f'SUCCESS: solved instance {num}!')
            df.loc[num,'solved']=1
            newly_solved.append(num)
        if int(opt)!=prev_opt:
            print(f'WARNING: found other solution for {num} ({opt} instead of {prev_opt})')
            df.loc[num,'opt'] = opt
        print(f'lower bound {low}')
        print(f'opt         {opt}')
        print(f's. space    {spa}')
        print(f'time        {tme} ms')
        print(f'diff space  {int(spa) - prev_space}')
        print(f'diff time   {int(tme) - prev_time} ms')
        df.loc[num,'space'] = int(spa)
        df.loc[num,'time'] = int(tme)
        df.loc[num,'low'] = int(low)
    else:
        if prev_sol:
            print(f'WARNING: did not solve instance {num}!')
            df.loc[num,'solved'] = 0
            df.loc[num,'space'] = 0
            df.loc[num,'time'] = 0
            newly_lost.append(num)
    print(flush=True)

print('num solved', df.solved.sum())
print('new solved', newly_solved)
print('new lost', newly_lost)

df.to_csv(f'results/{date}.csv')

