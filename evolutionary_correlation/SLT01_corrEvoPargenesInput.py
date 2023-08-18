# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 17:01:36 2021

@author: Administrator
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import os
import sys
import subprocess

def trimids(fasta):
    infile = r'./' + indir + r'/' + fasta
    outfile = r'./' + tempdir + r'/' + fasta
    with open(infile,'r') as source:
        with open(outfile,'w') as target:
            for line in source:
                if line[0] == '>':
                    target.write(line.split(' ')[0] + '\n')
                else:
                    target.write(line)

def unique_faa(fasta):
    infile = r'./' + tempdir + r'/' + fasta
    outfile = r'./' + outdir + r'/' + fasta
    subprocess.run('seqkit rmdup '
                   + infile
                   + ' -o '
                   + outfile, shell = True)

indir = sys.argv[1]
outdir = sys.argv[2]
tempdir = 'temp'
os.mkdir(tempdir)
files = [f for f in os.listdir(indir) if f[-3:] == 'faa']
if len(sys.argv) > 3:
    threads = sys.argv[3]
else:
    threads = cpu_count() - 1

if __name__ == '__main__':
    with Pool(threads) as p:
        p.map(trimids, files)

files = [f for f in os.listdir(tempdir) if f[-3:] == 'faa']

if __name__ == '__main__':
    with Pool(threads) as p:
        p.map(unique_faa, files)

os.rmdir(tempdir)