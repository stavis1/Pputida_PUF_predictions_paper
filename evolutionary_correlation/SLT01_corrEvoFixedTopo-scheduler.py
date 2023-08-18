# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 12:18:14 2021

@author: Administrator
"""
import sys
import os
from multiprocessing import cpu_count
from multiprocessing import Pool
import subprocess
from time import sleep
from Bio import Phylo
import re
from copy import deepcopy

def raxmlng(gene):
    out_prefix = out_dir + gene + '/' + gene
    raxmlng_logfile = mlsearch_run + gene + '/' + gene + '.raxml.log'
    with open(raxmlng_logfile, 'r') as logfile:
        loglines = logfile.readlines()
    run_cmd = loglines[next(i for i,l in enumerate(loglines) if l.startswith('RAxML-NG was called at')) + 2].strip()
    run_cmd = re.sub(r'--prefix [\w/]+', r'--prefix ' + out_prefix, run_cmd)
    run_cmd = re.sub(r'--seed 1 ', '', run_cmd)
    run_cmd = run_cmd + '--tree-constraint ' + out_prefix + '_constriant.newick'
    process = subprocess.Popen(run_cmd, shell = True)
    return process

def trim_tree(gene):
    gene_tree_newick = mlsearch_run + gene + '/' + gene + '.raxml.bestTree'
    gene_tree = Phylo.read(gene_tree_newick, 'newick')
    gene_taxa = [t.name for t in gene_tree.get_terminals()]
    constraint_tree = deepcopy(species_tree)
    species_taxa = [t.name for t in constraint_tree.get_terminals()]
    [constraint_tree.prune(t) for t in species_taxa if t not in gene_taxa]
    os.mkdir(out_dir + gene)
    Phylo.write(constraint_tree, out_dir + gene + '/' + gene + '_constriant.newick' , 'newick')
    return

pgene_dir = sys.argv[1]
out_dir = sys.argv[2]
if len(sys.argv) > 3:
    cores = int(sys.argv[3])
else:
    cores = cpu_count()

mlsearch_run = pgene_dir + r'mlsearch_run/results/'
species_tree = Phylo.read(pgene_dir + 'astral_run/output_species_tree.newick', 'newick')
queue = [f for f in os.listdir(mlsearch_run) if f[-3:] == 'faa']

with Pool(cores) as p:
    constriant_trees = p.map(trim_tree, queue)

running = []
while queue:
    for proc in running:
        if proc.poll() is not None:
            running.remove(proc)
    if len(running) < cores:
        running.append(raxmlng(queue[-1]))
        queue.pop(-1)
    else:
        sleep(1)
for proc in running:
    proc.wait()




