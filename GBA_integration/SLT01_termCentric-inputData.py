# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:20:11 2022

@author: Administrator
"""

import pandas as pd
import sys
import os
import dill
from collections import defaultdict

data_file = sys.argv[1]
protein_file = sys.argv[2]
term_file = sys.argv[3]
test_file = sys.argv[4]
out_dir = os.path.abspath(sys.argv[5]) + '/'

pair_data = pd.read_csv(data_file,sep = '\t')

with open(term_file,'rb') as dill_handle:
    terms = dill.load(dill_handle)
with open(protein_file,'rb') as dill_handle:
    proteins = dill.load(dill_handle)
unannotated = set([p.pp_id for p in proteins if p.unannotated])
prot_dict = {p.pp_id:p for p in proteins}
term_depth = {t.id:t.depth_long for t in terms}
good_go = set([t.id for t in terms])
for p in proteins:
    p.annotation = set([t for t in p.annotation if t in good_go])


with open(test_file,'r') as test_handle:
    test_set = set(test_handle.read().split('\n'))
print('read', flush=True)

def num_nonzero(x):
    return len([i for i in x if i != 0])

metrics = {'bitscore':(sum,max),'coEx':(sum,min),'coexpression_transferred':(max,sum),'combined_score':(sum,max),'cooccurence':(max,sum),
           'database_transferred':(sum,max),'database':(sum,max),'exJaccard':(sum,),'experiments_transferred':(sum,max),'fusion':(max,sum),
           'homology':(sum,max),'jaccard':(sum,min),'max_tm':(sum,max),'neighborhood_transferred':(sum,max),'neighborhood':(sum,max),'rcw':(sum,max),
           'operon':(num_nonzero,),'shared_taxa':(sum,),'textmining_transferred':(sum,),'textmining':(sum,max)}
metric_ids = list(metrics.keys())

accumulator = defaultdict(lambda : [[] for _ in range(len(metric_ids))])

def accumulate_data(row):
    for term in prot_dict[row[-1]].annotation:
        _=[accumulator[(row[-2],term)][i].append(row[i]) for i in range(len(metric_ids))]
    for term in prot_dict[row[-2]].annotation:
        _=[accumulator[(row[-1],term)][i].append(row[i]) for i in range(len(metric_ids))]
    return None

_ = [accumulate_data(row) for row in zip(*[pair_data[c] for c in metric_ids + ['protein1','protein2']])]
print('accumulated', flush=True)

def calc_metrics(key):
    features = [f(accumulator[key][i]) for i,metric in enumerate(metric_ids) for f in metrics[metric]]
    features.append(len(accumulator[key][0]))
    return features

term_data = pd.DataFrame({';'.join(k):calc_metrics(k) for k in accumulator.keys()}, index = [m+'-'+f.__name__ for m in metric_ids for f in metrics[m]] + ['n_hits'])
term_data = term_data.T
print('calculated', flush=True)


def label(protein, term):
    if prot_dict[protein].unannotated:
        return -1
    elif term_depth[term] > max([term_depth[t] for t in prot_dict[protein].annotation]):
        return -1
    elif term in prot_dict[protein].annotation:
        return 1
    else:
        return 0

def partition(protein):
    if prot_dict[protein].unannotated:
        return 'Q'
    elif protein in test_set:
        return 'Te'
    else:
        return 'Tr'

term_data['protein'] = [i.split(';')[0] for i in term_data.index]
term_data['term'] = [i.split(';')[1] for i in term_data.index]
term_data = term_data[[t in good_go for t in term_data['term']]]

term_data['label'] = [label(p,t) for p,t in zip(term_data['protein'],term_data['term'])]
term_data['partition'] = [partition(p) for p in term_data['protein']]

term_data.to_csv(f'{out_dir}termCentricFeatures.tsv',sep = '\t', index = False)
print('saved', flush=True)








