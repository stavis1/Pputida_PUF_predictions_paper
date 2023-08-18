# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 13:28:11 2022

@author: Administrator
"""

import sys
import os
import multiprocessing
import re

def good_spectrum(spectrum):
    name = re.search(r'Name: ([A-Z]+)', spectrum).group(1)
    if name in peps:
        return None
    elif any(re.search(name,seq) for seq in seqs):
        return None
    else:
        spectrum = re.sub(r'(Name: [^/]+/\d+_\d+)[^\n]+(\n)',r'\1\2',spectrum)
        return spectrum

def digest(seq):
    return [s for s in re.findall(r'([^KR]+[KR]|\Z)',seq) if len(s) > 6]

in_msp = os.path.abspath(sys.argv[1])
out_msp = os.path.abspath(sys.argv[2])
fasta = os.path.abspath(sys.argv[3])
cores = int(sys.argv[4]) if len(sys.argv) > 4 else multiprocessing.cpu_count()

peps = []
seqs = []
with open('Pseudomonas_putida_KT2440_110.faa','r') as faa:
    for line in faa:
        if not line.startswith('>'):
            seqs.append(line.strip())
            peps.append(digest(line.strip()))
peps = set([peptide for sublist in peps for peptide in sublist])

with open(in_msp,'r') as file:
    spectra = ['Name:' + s for s in file.read().split('Name:') if len(s) > 1]

if __name__ == '__main__':
    with multiprocessing.Pool(cores) as p:
        spectra = [s for s in p.map(good_spectrum, spectra) if s is not None]

with open(out_msp,'w') as file:
    for spectrum in spectra:
        file.write(spectrum)


