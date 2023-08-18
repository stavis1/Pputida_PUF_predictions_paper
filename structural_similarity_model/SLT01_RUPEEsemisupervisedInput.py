# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:17:34 2022

@author: Administrator
"""
import os
import pandas as pd
import dill
import networkx as nx
import numpy as np
from collections import defaultdict
import re
import sys
# from collections import Counter

in_dir = os.path.abspath(sys.argv[1]) + '/'
out_dir = os.path.abspath(sys.argv[2]) + '/'
puf_dir = os.path.abspath(sys.argv[3]) + '/'
ann_dir = os.path.abspath(sys.argv[4]) + '/'

#set up background information about Pseudomonas putida proteins and GO terms
with open(f'{in_dir}go_g.dill','rb') as pickle:
    go_g = dill.load(pickle)
with open(f'{in_dir}terms.dill','rb') as pickle:
    terms = dill.load(pickle)
with open(f'{in_dir}proteins.dill','rb') as pickle:
    proteins = dill.load(pickle)
prot_dict = {p.pp_id:p for p in proteins}
pp_uni = {p.pp_id:p.uni_id for p in proteins}
go_dict = {t.id:t for t in terms}

######
rng = np.random.default_rng(0)
tm_filter = 0.3
######

def all_parents(term):
    #get all terms reachable from query terms
    return list(nx.descendants(go_g,term))

def full_ann(termset):
    #for a set of terms include all terms that are also necessarily annotated due to the structure of the GO graph
    terms = []
    for term in termset:
        if term in valid_terms:
            terms.extend(all_parents(term.strip()))
    return set(terms)

def phylo_dist(lineage):
    #crude quantification of phylogenetic distance from Pseudomonas putida
    if 'Pseudomonas' in lineage:
        return 1
    elif 'Pseudomonadaceae' in lineage:
        return 2
    elif 'Pseudomonadales' in lineage:
        return 3
    elif 'Gammaproteobacteria' in lineage:
        return 4
    elif 'Proteobacteria' in lineage:
        return 5
    elif 'Bacteria' in lineage:
        return 6
    else:
        return 7

def tfu(pp,term,partition):
    #true/false/unannotated
    #sets up indicator for whether query term is actually annotated to protein
    if prot_dict[pp].unannotated:
        return -1
    elif term in prot_dict[pp].annotation:
        return 1
    else:
        return 0

def tetrq(pp):
    #test/train/query
    #partitions the annotated dataset into a 90/10 test train split
    if prot_dict[pp].unannotated:
        return 'Q'
    elif rng.random() > 0.1:
        return 'Tr'
    else:
        return 'Te'


puf_results = [f for f in os.listdir(puf_dir) if f.endswith('.txt')]
ann_results = [f for f in os.listdir(ann_dir) if f.endswith('.txt')]

#data structure to hold data on putative annotations
#[[tm-scores],[rmsds],[sequence_identitites],[sequence_gaps],[eukaryote],[pylogenetic_distance],n_hits,Te/Tr/Q,T/F/U]
term_data = defaultdict(lambda: [[],[],[],[],[],[],0,None,None])

#set up protein level data for hit proteins
uni_pattern = r'(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})|(?:[OPQ][0-9][A-Z0-9]{3}[0-9])'
valid_terms = set(t.id for t in terms)
all_uni_data = pd.read_csv(f'{in_dir}uniprot_hits.tsv', sep = '\t')
all_uni_data['uniprot'] = [re.search(uni_pattern,u).group() for u in all_uni_data['uniprot']]
uni_go = {u:full_ann(t.split(';')) if type(t) == str else None for u,t in zip(all_uni_data['uniprot'],all_uni_data['GO'])}
uni_euk = {u:'Eukaryota' in t for u,t in zip(all_uni_data['uniprot'],all_uni_data['lineage'])}
uni_phylo = {u:phylo_dist(t.split(';')) for u,t in zip(all_uni_data['uniprot'],all_uni_data['lineage'])}

#collect data from pairwise alignments, phylogeny, and set up indicators for train/test/query and true/false/unannotated
with open(f'{in_dir}pairwise_results.tsv','r') as alignment_file:
    for line in alignment_file:
        values = line.strip().split('\t')
        pp = values[0]
        uni = re.search(uni_pattern,values[1]).group()
        percent_id = float(values[2])
        percent_gap = float(values[3])
        donor_terms = uni_go[uni]
        if not donor_terms:
            continue
        for term in donor_terms:
            term_data[(pp,term)][2].append(percent_id)
            term_data[(pp,term)][3].append(1-percent_gap)
            term_data[(pp,term)][4].append(uni_euk[uni])
            term_data[(pp,term)][5].append(uni_phylo[uni])
            term_data[(pp,term)][6] += 1
            if term_data[(pp,term)][7] is None:
                term_data[(pp,term)][7] = tetrq(pp)
            if term_data[(pp,term)][8] is None:
                term_data[(pp,term)][8] = tfu(pp,term,term_data[(pp,term)][8])

#set up mapping between pdb chain IDs and uniprot accessions
hit_map = pd.read_csv(f'{in_dir}chain_map.tsv', sep = '\t')
hit_map['uniprot'] = [re.search(uni_pattern,u).group() for u in hit_map['uniprot']]
hit_dict = {c:u for c,u in zip(hit_map['pdb'],hit_map['uniprot'])}
good_chains = set(hit_dict.keys())


#RUPEE output for protein of unknown function searches
for result in puf_results:
    pp = re.search(r'PP_\d{4}', result).group()
    pp_ipro = []
    hits = pd.read_csv(puf_dir + result)
    hits['Chain Id']  = [p.upper()[:4] + ':' + p.upper()[4:] for p in hits['Chain Id']]
    hits = hits[[c in good_chains for c in hits['Chain Id']]]
    hits = hits[hits['TM-Score'] > tm_filter]
    for chain, tm, rmsd in zip(hits['Chain Id'], hits['TM-Score'], hits['RMSD']):
        donor_terms = uni_go[hit_dict[chain]]
        if not donor_terms:
            continue
        for term in donor_terms:
            term_data[(pp,term)][0].append(tm)
            term_data[(pp,term)][1].append(rmsd)
    
#RUPEE output for annotated protein searches
for result in ann_results:
    pp = re.search(r'PP_\d{4}', result).group()
    hits = pd.read_csv(ann_dir + result)
    hits['Chain Id']  = [p.upper()[:4] + ':' + p.upper()[4:] for p in hits['Chain Id']]
    hits = hits[[c in good_chains for c in hits['Chain Id']]]
    hits = hits[hits['TM-Score'] > tm_filter]
    for chain, tm, rmsd in zip(hits['Chain Id'], hits['TM-Score'], hits['RMSD']):
        donor_terms = uni_go[hit_dict[chain]]
        if not donor_terms:
            continue
        for term in donor_terms:
            term_data[(pp,term)][0].append(tm)
            term_data[(pp,term)][1].append(rmsd)

#from collected data generate predictors and write to file
header = ['tm_max','tm_sum','rmsd_min','rmsd_mean',
          'seqid_max','seqid_sum','notgap_max','notgap_sum',
          'eukaryote_any','eukaryote_all','phylo_min','phylo_mean','n_hits',
          'partition','label','protein','term']
with open(f'{out_dir}semiSupervisedInput-TMfiltered.tsv','w') as output_file:
    output_file.write('\t'.join(header) + '\n')
    for example in term_data.keys():
        if term_data[example][0]:
            pp = example[0]
            term = example[1]
            tm_max = max(term_data[example][0])
            tm_sum = sum(term_data[example][0])
            rmsd_min = min(term_data[example][1])
            rmsd_mean = np.mean(term_data[example][1])
            seqid_max = max(term_data[example][2])
            seqid_sum = sum(term_data[example][2])
            notgap_max = max(term_data[example][3])
            notgap_sum = sum(term_data[example][3])
            eukaryote_any = int(any(term_data[example][4]))
            eukaryote_all = int(all(term_data[example][4]))
            phylo_min = min(term_data[example][5])
            phylo_mean = np.mean(term_data[example][5])
            n_hits = term_data[example][6]
            partition = term_data[example][7]
            label = term_data[example][8]
            output_file.write('\t'.join([str(v) for v in [tm_max,tm_sum,rmsd_min,rmsd_mean,
                                                          seqid_max,seqid_sum,notgap_max,notgap_sum,
                                                          eukaryote_any,eukaryote_all,phylo_min,phylo_mean,n_hits,
                                                          partition,label,pp,term]]) + '\n')






