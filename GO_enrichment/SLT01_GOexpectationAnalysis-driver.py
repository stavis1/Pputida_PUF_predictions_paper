# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 16:15:30 2022

@author: Administrator
"""

from cmdstanpy import CmdStanModel
import pandas as pd
import os
import dill
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from collections import defaultdict
import networkx as nx

os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01')

with open('pickled_objects/proteins.dill','rb') as pickle:
    proteins = dill.load(pickle)
prot_dict = {p.pp_id:p for p in proteins}
with open('pickled_objects/terms.dill','rb') as pickle:
    terms = dill.load(pickle)
go_dict = {t.id:t for t in terms}
good_go = set(go_dict.keys())
with open('pickled_objects/go_g.dill','rb') as pickle:
    go_g = dill.load(pickle)

############ background data
#frequency of each GO term among PKFs
term_freqs = Counter([t for prot in proteins for t in prot.annotation if t in good_go])
#collect GO terms by their depth
depth_dict = defaultdict(lambda: [])
_=[depth_dict[go_dict[t].depth_long].append(t) for t in term_freqs.keys()] 
#frequency of GO term depths among PKFs
depth_freqs = {i:sum([term_freqs[t] for t in depth_dict[i]]) for i in range(max(depth_dict.keys()) + 1)}
#log odds of GO terms in PKFs
alpha_priors = {t:np.log(term_freqs[t]/(depth_freqs[go_dict[t].depth_long] - term_freqs[t])) for t in term_freqs.keys()}

############ PUF predictions
#fills in all missing parent terms for a set of GO annotations
def full_ann(termset):
    termset = [t for t in termset if t in good_go]
    all_terms = termset
    _=[all_terms.extend(nx.descendants(go_g,t)) for t in termset]
    return set(all_terms)

#guilt by association model data
gba_model = pd.read_csv('termCentric/predictions.tsv', sep = '\t')
gba_model = gba_model[[p == 'Q' for p in gba_model['partition']]]
gba_model = gba_model[gba_model['predicted']]

#protein structure similarity model data
structure_model = pd.read_csv('alphafold/predictions.tsv', sep = '\t')
structure_model = structure_model[[p == 'Q' for p in structure_model['partition']]]
structure_model = structure_model[structure_model['predicted']]

#collect model predictions
preds_accum = defaultdict(lambda : set([]))
_=[preds_accum[p].add(t) for p,t in zip(gba_model['protein'],gba_model['term'])]
_=[preds_accum[p].add(t) for p,t in zip(structure_model['protein'],structure_model['term'])]
preds = {p:full_ann(preds_accum[p]) for p in preds_accum.keys()}

#frequencies of GO terms in PUF predictions
pred_freqs = Counter([t for puf in preds.keys() for t in preds[puf]])
#collect predicted GO terms by their depth
pred_depth_dict = defaultdict(lambda: [])
_=[pred_depth_dict[go_dict[t].depth_long].append(t) for t in pred_freqs.keys()]
#frequency of GO term depths in PUF predictions
pred_depth_freqs = {i:sum([pred_freqs[t] for t in pred_depth_dict[i]]) for i in range(max(pred_depth_dict.keys()) + 1)}

#only estimate odds ratios for terms with at least 6 predicted observations
estimable = [p for p in pred_freqs.keys() if pred_freqs[p] > 3]

#compile executable stan model from file
stan_file = 'scripts/SLT01_GOexpectationAnalysis.stan'
model = CmdStanModel(stan_file=stan_file)

#input to stan model
data = {'N':len(estimable),
        'N_PUFs':len([p for p in proteins if p.unannotated]),
        'N_PUF_A':[pred_freqs[t] for t in estimable], #number of PUFs annotated with term t
        'N_PUF':[pred_depth_freqs[go_dict[t].depth_long] for t in estimable], #number of annotations of depth of term t in PUFs
        'PKF_alpha':[alpha_priors[t] for t in estimable],
        'PKF_confidence':.5} #number of annotations of depth of term t in PKFs
fit = model.sample(data=data)
draws = fit.draws_pd()
summary = fit.summary()

# #collect MCMC samples of odds ratios into a data frame
# ORs = {t:draws['odds_ratio[{n}]'.format(n = i + 1)] for i,t in enumerate(estimable)}
# odds_ratios = pd.DataFrame(ORs)
# sorted_cols = [c[1] for c in sorted([(np.median(odds_ratios[t]),t) for t in odds_ratios.columns])]

expects = draws[[c for c in draws.columns if c.startswith('expected_pop')]]
expects.columns = estimable
sorted_cols = [c[1] for c in sorted([(np.median(expects[t]),t) for t in expects.columns])]

sorted_cols = sorted_cols[:-3]
#plot horizontal violin plots
fig, ax = plt.subplots()
violins = ax.violinplot(expects[sorted_cols],
                        showmeans = False,
                        showmedians = False,
                        showextrema = False,
                        widths = 0.8,
                        vert = False)
for violin in violins['bodies']:
    violin.set_facecolor('k')
    violin.set_alpha(.8)
ax.scatter([np.median(expects[t]) for t in sorted_cols],
           range(1,len(sorted_cols) + 1),
           s = 1,
           c = 'white')
for i,term in enumerate(sorted_cols):
    lower, upper = np.quantile(expects[term], [.1,.9])
    ax.plot([lower, upper],[i+1, i+1],'-w', linewidth = .5)
ax.scatter([pred_freqs[t] for t in sorted_cols],
           [i + 1 for i in range(len(sorted_cols))],
           s = 1, c = 'r', label = 'prediction counts')
ax.set_yticks(np.asarray(range(len(sorted_cols))) + 1,
              [go_dict[t].name for t in sorted_cols],
              fontsize = 5)
ax.set_xlim((0,250))
ax.set_ylim((0,len(sorted_cols) + 1))
ax.set_box_aspect(1)
ax.set_xlabel('expected number of PUFs with annotation', fontsize = 8)
ax.legend(loc = 'lower right', prop={'size':8})
# fig.savefig('presentations_and_reports/paper_figures/GOexpectations.png')



