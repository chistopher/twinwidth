#!/usr/bin/python3

import glob
import subprocess

# loop over all files in the directory
for filename in glob.glob('*'):
    print(filename)
    subprocess.run(['xz', '-dk', filename])
