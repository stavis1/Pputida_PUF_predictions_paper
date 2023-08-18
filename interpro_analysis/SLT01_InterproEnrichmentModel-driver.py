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
from scipy.stats import gaussian_kde
from matplotlib import cm

os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01')

with open('pickled_objects/proteins.dill','rb') as pickle:
    proteins = dill.load(pickle)
prot_dict = {p.pp_id:p for p in proteins}

#
pufs = set([p.pp_id for p in proteins if p.unannotated])
pkfs = set([p.pp_id for p in proteins if not p.unannotated])

signalp = pd.read_csv('jaccard/Pseudomonas_putida_KT2440_110.gff3', sep = '\t', header=None, skiprows=1)
signalp['sig'] = signalp[2] + signalp[8]
pkf_couter = Counter(signalp[[p in pkfs for p in signalp[0]]]['sig'])
puf_couter = Counter(signalp[[p in pufs for p in signalp[0]]]['sig'])

interpro = pd.read_csv('jaccard/features.csv')
pkf_couter.update(interpro[[p in pkfs for p in interpro['Locus Tag']]]['External Signature Accession'])
puf_couter.update(interpro[[p in pufs for p in interpro['Locus Tag']]]['External Signature Accession'])

termlist = list(set(puf_couter.keys()).intersection(pkf_couter.keys()))
cutoff = 4
termlist = [t for t in termlist if puf_couter[t] > cutoff  and pkf_couter[t] > cutoff]
annotated = set(signalp[0]).union(interpro['Locus Tag'])

#compile executable stan model from file
stan_file = 'scripts/SLT01_InterproEnrichmentModel.stan'
model = CmdStanModel(stan_file=stan_file)

#input to stan model
data = {'N':len(termlist),
        'h_mu_mu':0, #hyperprior mean for mean of alpha distribution
        'h_mu_sig':0.1, #hyperprior standard dev for mean of alpha distribution
        'h_sig':0.1, #hyperprior scale param for standard dev of alpha distribution
        'N_PUF_A':[puf_couter[t] for t in termlist], #number of PUFs annotated with term t
        'N_PUF':[len([p for p in annotated if prot_dict[p].unannotated])]*len(termlist), #number of annotations of depth of term t in PUFs
        'N_PKF_A':[pkf_couter[t] for t in termlist], #number of PKFs annotated with term t
        'N_PKF':[len([p for p in annotated if not prot_dict[p].unannotated])]*len(termlist)} #number of annotations of depth of term t in PKFs
fit = model.sample(data=data)
draws = fit.draws_pd()
summary = fit.summary()

#collect MCMC samples of odds ratios into a data frame
ORs = {t:draws['odds_ratio[{n}]'.format(n = i + 1)] for i,t in enumerate(termlist)}
odds_ratios = pd.DataFrame(ORs)
sorted_cols = [c[1] for c in sorted([(np.median(odds_ratios[t]),t) for t in odds_ratios.columns])]
label_names = {'SSF54909':'Dimeric alpha+beta barrel',
 'G3DSA:3.40.50.300':'P-loop containing nucleotide triphosphate hydrolases',
 'PS51257':'Prokaryotic membrane lipoprotein lipid attachment site profile',
 'PF13007':'Transposase C of IS166 homeodomain',
 'mobidb-lite':'consensus disorder prediction',
 'PF03050':'Transposase IS66 family',
 'G3DSA:3.30.70.100':'CATH Superfamily 3.30.70.100',
 'signal_peptide.':'signal peptide',
 'G3DSA:3.40.190.10':'Periplasmic binding protein-like II (Gene3D)',
 'SSF53850':'Periplasmic binding protein-like II (Superfamily 1.75)',
 'signal_peptideNote=TAT':'TAT signal peptide',
 'G3DSA:3.10.450.50':'CATH Superfamily 3.10.450.50',
 'PF13005':'zinc-finger binding domain of transposase IS66',
 'lipoprotein_signal_peptide.':'lipoprotein signal peptide',
 'PF13817':'IS66 C-terminal element',
 'Coil':'Coil'}
def scatter(scale):
    return (np.random.random()*(scale/2)) - (scale/4)    

def dense_color(values):
    values = np.log(np.asarray(values))
    density_obj = gaussian_kde(values)
    densities = density_obj.evaluate(values)
    low = min(densities)
    high = max(densities)
    return [cm.plasma(int(((val-low)/(high-low))*cm.plasma.N)) for val in densities]

#plot horizontal violin plots
fig, axes = plt.subplots(1,2, layout = 'constrained')
fig.set_figheight(6)
fig.set_figwidth(10)
ax = axes[0]
ax.grid(axis='y')
ax.set_axisbelow(True)
for i,term in enumerate(sorted_cols):
    ax.scatter(odds_ratios[term],
                [i+ 1 + scatter(1.75) for _ in range(odds_ratios.shape[0])],
                color = dense_color(odds_ratios[term]),
                s = 1, marker = '.', alpha = 0.1)
    lower, upper = np.quantile(ORs[term], [.1,.9])
    ax.plot([lower, upper],[i+1, i+1],'-k', linewidth = 1)
ax.scatter([np.median(ORs[t]) for t in sorted_cols],
           range(1,odds_ratios.shape[1] + 1),
           s = 3,
           c = 'k')
ax.set_xscale('log')
ax.plot([1,1],[0,odds_ratios.shape[1] + 1],'--r', linewidth = .5)
ax.set_yticks(np.asarray(range(len(sorted_cols))) + 1,
              [label_names[c] for c in sorted_cols],
              fontsize = 7)
ax.set_xticks([.1,1,10,100],['1/10','1/1','10/1','100/1'])
ax.set_xlabel('Posterior Odds Ratio')
ax.set_ylim((0,odds_ratios.shape[1] + 1))
ax.set_box_aspect(1)

ax = axes[1]
# ax.set_ylim((0,odds_ratios.shape[1] + 1))
ax.barh(y = range(len(sorted_cols)), 
        width = [puf_couter[t] for t in sorted_cols], 
        height = 0.4, align = 'edge', color = 'k',
        label = 'PUF counts')
ax.barh(y = range(len(sorted_cols)), 
        width = [pkf_couter[t] for t in sorted_cols], 
        height = -0.4, align = 'edge', color = 'r',
        label = 'PKF counts')
ax.set_xscale('log')
ax.set_ylim((-1,odds_ratios.shape[1]))
ax.set_box_aspect(1)
ax.set_yticks([])
ax.legend()
ax.set_xlabel('Counts')
fig.suptitle('InterProScan Element Enrichments', y = .85)

# plt.close('all')


