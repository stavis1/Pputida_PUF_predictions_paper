# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:33:02 2022

@author: Administrator
"""

import dill
import os
from sortedcontainers import SortedSet
from collections import defaultdict
import pandas as pd

os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01')

with open('pickled_objects/proteins.dill','rb') as pickle:
    proteins = dill.load(pickle)
with open('pickled_objects/terms.dill','rb') as pickle:
    terms = dill.load(pickle)

go_dict = {t.id:t for t in terms}
pp_uni = {p.pp_id:p.uni_id for p in proteins}
good_terms = set([t.id for t in terms])

############################## combine edgelists for GBA model ###########################

annotation_long = {p.pp_id:SortedSet([(go_dict[g].depth_long,g) for g in p.annotation if g in good_terms], key =lambda x: x[0]) for p in proteins}

def to_pp(name):
    name = name.strip()
    if name.startswith('PP_') and name in annotation_long.keys():
        return name
    elif name.startswith('PP_'):
        uni = pp_uni[name]
        return next(p.pp_id for p in proteins if uni == p.uni_id)
    else:
        return next(p.pp_id for p in proteins if name == p.uni_id)

def sim_metrics(pair):
    pair = [p.strip() for p in pair]
    shared_long = annotation_long[pair[0]].intersection(annotation_long[pair[1]])
    if len(shared_long) == 0:
        return 0
    else:
        shared_long = annotation_long[pair[0]].intersection(annotation_long[pair[1]])
        return shared_long[-1][0]

#cols holds all of the predictor names along with their null value
cols = {'coEx':0,'exJaccard':0,'shared_taxa':0,'rcw':100,'bitscore':0,'jaccard':1,
        'neighborhood':0,'neighborhood_transferred':0,'fusion':0,'cooccurence':0,'homology':0,
        'coexpression_transferred':0,'experiments_transferred':0,'operon':0,'rmsd':11,'max_tm':0,
        'database':0,'database_transferred':0,'textmining':0,'textmining_transferred':0,'combined_score':0}
collist = list(cols.keys())
pairlist = defaultdict(lambda: [cols[k] for k in collist])

for edgelist in [l for l in os.listdir('edgelists')]:
    with open('edgelists/' + edgelist,'r') as file:
        header = file.readline().split('\t')
        colmap = {}
        for i,col in enumerate(header[2:]):
            colmap[i] = collist.index(col.strip())
        for line in file:
            entries = line.strip().split('\t')
            if entries[0] != entries[1]:
                pair = frozenset([to_pp(entries[0]),to_pp(entries[1])])
                for i,entry in enumerate(entries[2:]):
                    pairlist[pair][colmap[i]] = entry

for pair in pairlist.keys():
    pairlist[pair].append(sim_metrics(pair))

sim_data = {c:[pairlist[p][i] for p in pairlist.keys()] for i,c in enumerate(collist)}
sim_data['depth_long'] = [pairlist[p][-1] for p in pairlist.keys()]
sim_data['protein1'] = [list(p)[0] for p in pairlist.keys()]
sim_data['protein2'] = [list(p)[1] for p in pairlist.keys()]
sim_dataframe = pd.DataFrame(sim_data)

del sim_data, pairlist

sim_dataframe.to_csv('model_data/SLT01_combinedEdgelist.tsv', sep = '\t', index = False)



