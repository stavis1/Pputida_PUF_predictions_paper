# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 12:11:49 2021

@author: Administrator
"""

import os
import sys
from multiprocessing import Pool
from multiprocessing import cpu_count
import subprocess

def mf_align(fasta):
    infile = r'./' + indist + r'/' + fasta
    outfile = r'./' + outdist + r'/' + fasta
    if os.path.isfile(outfile):
        return
    mafft_process = subprocess.run(['/usr/bin/mafft',
                                    '--auto',
                                    infile,
                                    '>',
                                    outfile])
    return

indist = sys.argv[1]
outdist = sys.argv[2]
orthogroups = [f for f in os.listdir(indist) if f[-3:] == 'faa']
if len(sys.argv) > 3:
    threads = int(sys.argv[3])
else:
    threads = cpu_count() - 1
if __name__ == '__main__':
    with Pool(threads) as p:
        p.map(mf_align, orthogroups)