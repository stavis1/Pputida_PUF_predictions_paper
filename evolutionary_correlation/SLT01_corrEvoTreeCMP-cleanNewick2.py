# -*- coding: utf-8 -*-
"""
Created on Thu May  6 13:53:54 2021

@author: Administrator
"""

import sys
import os
import re
from multiprocessing import Pool
from multiprocessing import current_process
from multiprocessing import cpu_count
from Bio import Phylo
from io import StringIO

def dict_data_setup(fasta):
    outfile = newick_dir + 'temp' + str(current_process().pid) + '.txt'
    idlist = []
    with open(fasta_dir + fasta,'r') as faa:
        for line in faa:
            if line[0] == '>':
                idlist.append(line.split(' ')[0][1:] + ',' + fasta[:-4] + '\n')
    if os.path.exists(outfile):
        with open(outfile,'a') as out:
            out.writelines(idlist)
    else:
        with open(outfile,'w') as out:
            out.writelines(idlist)

def clean_newick(protein):
    newick_file = results_dir + protein + '/' + protein + '.raxml.bestTree'
    clean_newick_file = newick_dir + protein + '.newick'
    with open(newick_file) as newick:
        tree = newick.readline()
    orthologs = re.findall(gene_pattern, tree)#gene names from the newick string
    for ortholog in orthologs:
        tree = re.sub(ortholog, prot_strain_dict[ortholog], tree)#replaces gene names with strain names for TreeCMP and ASTRAL
    tree = Phylo.read(StringIO(tree), 'newick')#converts tree from string to Phylo tree object
    taxa = [t.name for t in tree.get_terminals()]
    dup_taxa = list(set([t for t in taxa if taxa.count(t) > 1]))
    pp_gene = tree.find_any(name = 'Pseudomonas_putida_KT2440_110')
    for taxon in dup_taxa:
        genes = list(tree.find_elements(terminal = True, name = taxon))#gets clade objects for each duplicate taxon
        min_dist = min([tree.distance(pp_gene, g) for g in genes])#finds the minimun distance between pputida and a duplicate taxon
        gene_to_save = [g for g in genes if tree.distance(pp_gene, g) == min_dist][0]
        [tree.prune(g) for g in genes if g != gene_to_save]#prune operation automatically removes internal nodes with one descendent
    Phylo.write(tree, clean_newick_file, 'newick')
    return



fasta_dir = sys.argv[1]
results_dir = sys.argv[2]
newick_dir = sys.argv[3]
if len(sys.argv) > 4:
    threads = sys.argv[4]
else:
    threads = cpu_count()

gene_pattern = re.compile(r'[a-zA-Z]\w+')


fastas = os.listdir(fasta_dir)
if __name__ == '__main__':
    with Pool(threads) as p:
        p.map(dict_data_setup, fastas)

tempfiles = os.listdir(newick_dir)

prot_strain_dict = {}
for tempfile in tempfiles:
    with open(newick_dir + tempfile, 'r') as file:
        for line in file:
            entry = line.split(',')
            prot_strain_dict[entry[0]] = entry[1]

for tempfile in tempfiles:
    os.unlink(newick_dir + tempfile)

proteins = os.listdir(results_dir)
if __name__ == '__main__':
    with Pool(threads) as p:
        p.map(clean_newick, proteins)







