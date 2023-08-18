# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:36:59 2022

@author: Administrator
"""

import os
import sys
import subprocess
from multiprocessing import Pool
import re

jobs_file = os.path.abspath(sys.argv[1])
uni_fasta = os.path.abspath(sys.argv[2])
pp_fasta = os.path.abspath(sys.argv[3])
out_file = os.path.abspath(sys.argv[4])
if len(sys.argv) > 5:
    cores = int(sys.argv[5])
else:
    cores = os.cpu_count()

def run_nw(pp_id,uni_id):
    #run NWalign through shell and wait for it to finish
    results = subprocess.run(' '.join(['NWalign',
                                      pp_seq_dict[pp_id],
                                      uni_seq_dict[uni_id],
                                      '3']),
                            shell = True,
                            capture_output = True,
                            text = True)
    
    #parse out the percentage identical aa and percentage gaps in the alignment
    percent_id = re.search(r'Sequence identity:\s+(\d\.\d+)',results.stdout).group(1)
    percent_gap = results.stdout.split('\n')[7].count('-')/len(results.stdout.split('\n')[7])
    return '\t'.join([str(v) for v in [pp_id,uni_id,percent_id,percent_gap]])

#get sequence information for Pseudomonas putida proteins
with open(pp_fasta,'r') as fasta:
    entries = fasta.read().split('\n>')
pp_seq_dict = {re.search(r'PP_\d{4}',e).group():''.join(e.split('\n')[1:]) for e in entries}

#get sequence information for proteins identified by RUPEE
uni_pattern = r'(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})|(?:[OPQ][0-9][A-Z0-9]{3}[0-9])'
with open(uni_fasta,'r') as fasta:
    entries = fasta.read().split('\n>')
uni_seq_dict = {re.search(uni_pattern,e).group():''.join(e.split('\n')[1:]) for e in entries}

#each job is a P. putida protein and a uniprot hit identified by RUPEE
with open(jobs_file,'r') as all_jobs:
    jobs = [j.split('\t') for j in all_jobs.read().split('\n') if j]

#ensure all proteins have sequence information
good_pp_ids = set(pp_seq_dict.keys())
assert all(j[0] in good_pp_ids for j in jobs)
good_uni_ids = set(uni_seq_dict.keys())
jobs = [j if j[1] in good_uni_ids else (j[0],re.search(uni_pattern,j[1]).group()) for j in jobs]
assert all(j[1] in good_uni_ids for j in jobs)

#run NWalign in parallel on all jobs
with Pool(cores) as p:
    results = p.starmap(run_nw,jobs)

#save results
with open(out_file,'w') as out:
    out.write('\n'.join(results))



