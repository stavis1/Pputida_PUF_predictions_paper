# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:08:47 2022

@author: Administrator
"""

import dill
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import re
from scipy.stats import gaussian_kde
import venn
# import itertools
import networkx as nx
from sortedcontainers import SortedList
from collections import Counter
from scipy.integrate import trapz

os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01_PUFs')

with open('pickled_objects/proteins.dill','rb') as pickle:
    proteins = dill.load(pickle)
prot_dict = {p.pp_id:p for p in proteins}
uni_pp = {p.uni_id:p.pp_id for p in proteins}
good_uni = set(uni_pp.keys())
with open('pickled_objects/terms.dill','rb') as pickle:
    terms = dill.load(pickle)
go_dict = {t.id:t for t in terms}
good_go = set(go_dict.keys())

#######Figure 1
#this one was made in inkscape

#######Figure 2

uni_anns = []
with open('baseline_information/names_seq_input.csv','r') as uni_data:
    for line in uni_data:
        loci = re.findall(r'PP_\d{4}',line)
        gos = re.findall(r'GO:\d{7}',line)
        uni_anns.extend([(l,g) for l in loci for g in gos])

bad_count = 0
biocyc_anns = []
with open('baseline_information/proteins.dat','r') as biocyc_data:
    bc_entries = biocyc_data.read().split('//')
for entry in bc_entries:
    try:
        gene = re.search(r'UNIPROT "([^"]+)"',entry).group(1)
        if gene in good_uni:
            locus = uni_pp[gene]
        for go in re.findall(r'GO:\d{7}', entry):
            biocyc_anns.append((locus,go))
    except:
        bad_count += 1

pseudDB = pd.read_csv('baseline_information/gene_ontology_csv.csv')
pseudDB_anns = [(l,g) for l,g in zip(pseudDB['Locus Tag'],pseudDB['Accession'])]

netgo_anns = []
netgo_files = [f'baseline_information/{f}' for f in os.listdir('baseline_information') if f.startswith('result_')]
for result_file in netgo_files:
    netgo = pd.read_csv(result_file, sep = '\t', header = None)
    netgo_anns.extend([(uni_pp[l],g) for l,g,p in zip(netgo[0],netgo[1],netgo[2]) if (float(p) > 0.9) and l in good_uni])

venn.venn({'Uniprot':set(uni_anns),
           'Biocyc':set(biocyc_anns),
           'Pseudomonas Database':set(pseudDB_anns),
           'NetGO2.0':set(netgo_anns)})
plt.savefig('presentations_and_reports/paper_figures/GOvenn.svg')


depths = [max([go_dict[t].depth_long for t in p.annotation if t in good_go]) if not p.unannotated else 0 for p in proteins]
depths_counter = Counter(depths)
plt.bar(range(max(depths) + 1), 
        [depths_counter[i] for i in range(max(depths) + 1)],
        color = 'k')
plt.ylabel('Number of Proteins')
plt.xlabel('Depth of Deepest GO Annotation')
plt.savefig('presentations_and_reports/paper_figures/GOdepths.svg')

#######Figure 3
def color(values):
    low = min(values)
    high = max(values)
    return [cm.plasma(int(((val-low)/(high-low))*cm.plasma.N)) for val in values]

def get_sm(vals):
    colormap = plt.cm.get_cmap('plasma')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin = min(vals), vmax = max(vals))
    return sm

#read in the output of the functional similarity prediction model and filter for expected similar proteins
similarity_cutoff = 6
prot_sim = pd.read_csv('termCentric/SLT01_proteinSimilarityModelResults.tsv',sep = '\t')
prot_sim = prot_sim[prot_sim['predicted_depth'] > similarity_cutoff]

#construct an unweighted, undirected networkx graph and include unconnected nodes 
sim_g = nx.Graph(zip(prot_sim['protein1'],prot_sim['protein2']))
sim_g.add_nodes_from([p.pp_id for p in proteins])


prots = set(sim_g.nodes())
pufs = set([p.pp_id for p in proteins if p.unannotated])

def calc_modularity(G, partition, gamma):
    return nx.algorithms.community.quality.modularity(G,partition,resolution = gamma)

#this procedure comes from https://doi.org/10.1103/PhysRevE.94.052315
#it's a way to estimate the optimal resolution parameter for modularity based on the relationship
#between modularity and the degree-corrected planted partition model of communities in networks.
#greedy_modularity_communities is a little slow so this may take a couple minutes
def calc_gamma(G, partition):
    e_total = len(G.edges())
    e_within = 0
    o_in_divisor = 0
    for part in partition:
        sub_g = nx.subgraph(G, part)
        e_within += len(sub_g.edges())
        o_in_divisor += (sum([d[1] for d in nx.degree(G,part)])**2)/(2*e_total)
    e_between = e_total - e_within
    omega_in = (2*e_within)/o_in_divisor
    omega_out = (2*e_between)/((2*e_total) - o_in_divisor)
    return (omega_in - omega_out)/(np.log(omega_in) - np.log(omega_out))
old_gamma = []
gamma = 1
flag = True
for _ in range(15):
    old_gamma.append(gamma)
    partition = nx.algorithms.community.greedy_modularity_communities(sim_g, resolution = gamma)
    gamma = calc_gamma(sim_g, partition)
    if abs(old_gamma[-1] - gamma) < 0.001:
        old_gamma.append(gamma)
        break

#the modularity of PUFs using the optimal resolution parameter
modularity = calc_modularity(sim_g, [pufs, prots.difference(pufs)], gamma)

#random partitions are drawn and their modularities are measured
prot_arr = np.asarray(list(prots))
null_modularities = []
null_degree_seqs = []
for _ in range(1000):
    part = np.random.choice(prot_arr, len(pufs), replace = False)
    null_modularities.append(calc_modularity(sim_g, [part, prots.difference(part)], gamma))
    null_degree_seqs.append(Counter([d[1] for d in list(nx.degree(sim_g, part))]))

#To contextualize the modularity scores we look at the modularity of sets of proteins all
#annotated with the same term. Because the balance of the partition makes a difference to the score
#we only look at terms that are annotated to a large number of proteins
term_counts = Counter([t for p in proteins for t in p.annotation])
dist_terms = [t.id for t in terms if term_counts[t.id] > 400]
dist_sizes = [term_counts[go_dict[t].id] for t in dist_terms]
dist_mods = []
for t in dist_terms:
    ann_prots = set([p for p in prots if t in prot_dict[p].annotation])
    dist_mods.append(calc_modularity(sim_g,[ann_prots, prots.difference(ann_prots)], gamma))

#Figure 3
fig,ax = plt.subplots()
ax.hist(null_modularities, bins = 50, color = 'r', alpha = 0.5, label = 'Random Modularity')
ax.hist(dist_mods, bins = 50, color = 'b', alpha = 0.5, label = 'Shared Function Modularity')
ylim = ax.get_ylim()
ax.plot([modularity]*2, ylim, '-k', linewidth = 0.5, label = 'PUF Modularity')
ax.set_ylim(ylim)
ax.legend(loc = 'upper center')
ax.set_xlabel('Modularity Score')
ax.set_ylabel('Partition Count')
fig.savefig('presentations_and_reports/paper_figures/Figure3.svg')


#Table S1
term_modularities = pd.DataFrame({'term':dist_terms,
                                  'term_name':[go_dict[t].name for t in dist_terms],
                                  'n_proteins':dist_sizes,
                                  'modularity':dist_mods})
term_modularities.sort_values(by = 'modularity', inplace = True)
term_modularities.to_csv('presentations_and_reports/paper_figures/TableS1.tsv',sep = '\t', index = False)


#######Figure 4
def scatter(scale):
    return (np.random.random()*(scale/2)) - (scale/4)    

def TPR(rates):
    divisor = rates[0] + rates[3]
    if divisor == 0:
        return 0
    else:
        return rates[0]/divisor

def FPR(rates):
    divisor = rates[1] + rates[2]
    if divisor == 0:
        return 0
    else:
        return rates[1]/divisor

def FDR(rates):
    tp,fp,tn,fn = rates
    divisor = tp + fp
    if divisor == 0:
        return 0
    else:
        return fp/divisor

def control_fdr(true_scores, false_scores):
    true_scores = SortedList(true_scores)
    false_scores = SortedList(false_scores)
    all_scores = true_scores + false_scores
    n_true = len(true_scores)
    n_false = len(false_scores)
    fdr_resids = []
    for score in all_scores:
        fn = true_scores.bisect_left(score)
        tp = n_true - fn
        tn = false_scores.bisect_left(score)
        fp = n_false - tn
        fdr_resids.append(abs(0.01 - FDR([tp,fp,tn,fn])))
    best_resid = min(fdr_resids)
    return next(c for c,r in zip(all_scores,fdr_resids) if r == best_resid)

def roc(true_scores, false_scores, fdr_cut, title, outfile):
    true_scores = SortedList(true_scores)
    false_scores = SortedList(false_scores)
    all_scores = true_scores + false_scores
    tprs = [0]
    fprs = [0]
    n_true = len(true_scores)
    n_false = len(false_scores)
    for score in reversed(all_scores):
        fn = true_scores.bisect_left(score)
        tp = n_true - fn
        tn = false_scores.bisect_left(score)
        fp = n_false - tn
        
        tprs.append(TPR([tp,fp,tn,fn]))
        fprs.append(FPR([tp,fp,tn,fn]))
    
    fn = true_scores.bisect_left(fdr_cut)
    tp = n_true - fn
    tn = false_scores.bisect_left(fdr_cut)
    fp = n_false - tn
    fdr_tpr = TPR([tp,fp,tn,fn])
    fdr_fpr = FPR([tp,fp,tn,fn])
    fdr = FDR([tp,fp,tn,fn])
    
    tprs.append(1)
    fprs.append(1)
    fig, ax = plt.subplots()
    ax.plot(fprs,tprs,color = 'k', linewidth = 1, clip_on=False)
    ax.set_box_aspect(1)
    ax.plot((0,1),(0,1),'--r', linewidth = .5)
    ax.set_ylabel('True Positive Rate')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylim((0,1.005))
    ax.set_xlim((-0.005,1))
    ax.scatter(fdr_fpr,fdr_tpr, c = 'r', s = 10, marker = '.', zorder=10)
    ax.text(fdr_fpr + 0.02, fdr_tpr -0.008, f'Realized FDR = {round(fdr,3)}', fontsize = 8)
    ax.set_title(title)
    fig.savefig(outfile)
    
    print(f'AUC: {trapz(tprs,fprs)}')
    print(f'recall at controlled FDR: {fdr_tpr}')
    return


gba_preds = pd.read_csv('termCentric/predictions.tsv', sep = '\t')
train_preds = gba_preds[gba_preds['partition'] == 'Tr']
train_true = train_preds[train_preds['label'] == 1]['prediction']
train_false = train_preds[train_preds['label'] == 0]['prediction']
fdr_cut = control_fdr(train_true, train_false)

puf_preds = gba_preds[gba_preds['partition'] == 'Q']
term_list = puf_preds[puf_preds['prediction'] > fdr_cut]['term']
gba_term_counts = Counter(term_list)
gba_depth_counts = Counter([go_dict[t].depth_long for t in term_list])

test_preds = gba_preds[gba_preds['partition'] == 'Te']
true_scores = test_preds[test_preds['label'] == 1]['prediction']
false_scores = test_preds[test_preds['label'] == 0]['prediction']
roc(true_scores, 
    false_scores, 
    fdr_cut, 
    'GBA Model Test Set',
    'presentations_and_reports/paper_figures/GBA_ROC.svg')


ss_preds = pd.read_csv('alphafold/predictions.tsv', sep = '\t')
train_preds = ss_preds[ss_preds['partition'] == 'Tr']
train_true = train_preds[train_preds['label'] == 1]['prediction']
train_false = train_preds[train_preds['label'] == 0]['prediction']
fdr_cut = control_fdr(train_true, train_false)

puf_preds = ss_preds[ss_preds['partition'] == 'Q']
term_list = puf_preds[puf_preds['prediction'] > fdr_cut]['term']
ss_term_counts = Counter(term_list)
ss_depth_counts = Counter([go_dict[t].depth_long for t in term_list])

test_preds = ss_preds[ss_preds['partition'] == 'Te']
true_scores = test_preds[test_preds['label'] == 1]['prediction']
false_scores = test_preds[test_preds['label'] == 0]['prediction']
roc(true_scores, 
    false_scores, 
    fdr_cut, 
    'Structural Similarity Model Test Set',
    'presentations_and_reports/paper_figures/SSim_ROC.svg')

deepest = max(max(gba_depth_counts.keys()),max(ss_depth_counts.keys()))
gba_bars = [gba_depth_counts[i] for i in range(deepest + 1)]
ss_bars = [ss_depth_counts[i] for i in range(deepest + 1)]

fig, ax = plt.subplots()
ax.bar(range(deepest + 1), gba_bars, color = 'r', label = 'GBA')
ax.bar(range(deepest + 1), ss_bars, color = 'k', bottom = gba_bars, label = 'Structural Similarity')
ax.set_xticks(range(deepest + 1))
ax.legend()
ax.set_ylabel('Number of Predicted Term Transfers')
ax.set_xlabel('Term Depth')
fig.savefig('presentations_and_reports/paper_figures/PredictedTermDepths.svg')

#######Figure 5

#######Figure 6

#######Figure S1
def scatter(scale):
    return (np.random.random()*(scale/2)) - (scale/4)    

def dense_color(values):
    values = np.log(np.asarray(values))
    density_obj = gaussian_kde(values)
    densities = density_obj.evaluate(values)
    low = min(densities)
    high = max(densities)
    return [cm.plasma(int(((val-low)/(high-low))*cm.plasma.N)) for val in densities]

def bootstrap_qq(x,y):
    qq_lines = []
    for _ in range(500):
        x_new = np.random.choice(x, size = len(x), replace = True)
        y_new = np.random.choice(y, size = len(y), replace = True)
        qq_lines.append((np.quantile(x_new,np.linspace(.01,.99,100)), 
                         np.quantile(y_new,np.linspace(.01,.99,100))))
    return qq_lines

uniprot_id_pattern = r'(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})|(?:[OPQ][0-9][A-Z0-9]{3}[0-9])'
with open('alphafold/SLT01_alphafoldPredsSummaryMetrics.txt','r') as summary:
    qual_list = summary.readlines()
qual_list = [q for q in qual_list if not re.search('G1JLS8',q)]
loci = [uni_pp[re.search(uniprot_id_pattern,q).group()] for q in qual_list]
mean_plddt = [float(q.split('\t')[1]) for q in qual_list]

pkf_quals = [m for m,l in zip(mean_plddt, loci) if not prot_dict[l].unannotated]
puf_quals = [m for m,l in zip(mean_plddt, loci) if prot_dict[l].unannotated]

fade = 0.5
fig,ax = plt.subplots()
ax.scatter([scatter(1) for _ in range(len(pkf_quals))], pkf_quals, 
            s = 1, color = 'k', marker = '.', alpha = fade)
ax.scatter([1 + scatter(1) for _ in range(len(puf_quals))], puf_quals, 
            s = 1, color = 'k', marker = '.', alpha = fade)
ax.set_ylabel('mean pLDDT')
ax.set_box_aspect(1)
ax.set_xticks((0,1),('PKFs','PUFs'))
fig.savefig('presentations_and_reports/paper_figures/alphafoldQualMetrics.svg')

low = min(np.quantile(pkf_quals,.01),np.quantile(puf_quals,.01))
high = max(np.quantile(pkf_quals,.99),np.quantile(puf_quals,.99))
fig,ax = plt.subplots()
ax.scatter(np.quantile(pkf_quals,np.linspace(.01,.99,100)),
           np.quantile(puf_quals,np.linspace(.01,.99,100)),
           s = 1, color = 'w', marker = '.', zorder=2)
ax.scatter(np.quantile(pkf_quals,.5),
           np.quantile(puf_quals,.5),
           s = 3, color = 'r', marker = '.', zorder=2)
ax.text(np.quantile(pkf_quals,.5) + 1, 
        np.quantile(puf_quals,.5) - 1, 
        'Median', fontsize = 6, zorder=2)
for i,line in enumerate(bootstrap_qq(pkf_quals, puf_quals)):
    ax.plot(line[0],line[1],'-k', linewidth = 0.2, alpha = 0.2, zorder=1)
ax.plot((low,high),(low,high),'--r',linewidth = 0.5)
ax.set_box_aspect(1)
ax.set_xlabel('PKF pLDDT')
ax.set_ylabel('PUF pLDDT')
ax.set_facecolor('lightgrey')
ax.set_title('Q-Q Plot of Mean pLDDT')
fig.savefig('presentations_and_reports/paper_figures/alphafoldQualMetrics-QQ.svg')


pkfs = [p for p in proteins if not p.unannotated]
pufs = [p for p in proteins if p.unannotated]
fade = 0.5
fig,ax = plt.subplots()
ax.scatter([scatter(1) for _ in range(len(pkfs))],
           [len(p.seq) for p in pkfs], 
           s = 1, color = 'k', marker = '.', alpha = fade)
ax.scatter([1 + scatter(1) for _ in range(len(pufs))],
           [len(p.seq) for p in pufs], 
           s = 1, color = 'k', marker = '.', alpha = fade)
ax.set_yscale('log')
ax.set_ylabel('Sequence Length')
ax.set_box_aspect(1)
ax.set_xticks((0,1),('PKFs','PUFs'))
fig.savefig('presentations_and_reports/paper_figures/seqLengths.svg')

pkf_lens = [len(p.seq) for p in pkfs]
puf_lens = [len(p.seq) for p in pufs]
short = min(np.quantile(pkf_lens,.01),np.quantile(puf_lens,.01))
long = max(np.quantile(pkf_lens,.99),np.quantile(puf_lens,.99))
fig,ax = plt.subplots()
ax.scatter(np.quantile(pkf_lens,np.linspace(.01,.99,100)),
           np.quantile(puf_lens,np.linspace(.01,.99,100)),
           s = 1, color = 'w', marker = '.', zorder=2)
ax.scatter(np.quantile(pkf_lens,.5),
           np.quantile(puf_lens,.5),
           s = 3, color = 'r', marker = '.', zorder=2)
ax.text(np.quantile(pkf_lens,.5) + 20, 
        np.quantile(puf_lens,.5) - 20, 
        'Median', fontsize = 6)
ax.plot((short,long),(short,long),'--r',linewidth = 0.5)
for i,line in enumerate(bootstrap_qq(pkf_lens, puf_lens)):
    ax.plot(line[0],line[1],'-k', linewidth = 0.2, alpha = 0.2, zorder=1)
ax.set_box_aspect(1)
ax.set_xlabel('PKF Lengths')
ax.set_ylabel('PUF Lengths')
ax.set_facecolor('lightgrey')
ax.set_title('Q-Q Plot of Sequence Lengths')
fig.savefig('presentations_and_reports/paper_figures/seqLengths-QQ.svg')


def n_orthologs(protein):
    with open(f'correlatedEvolution/orthogroupFaa/{protein}.faa','r') as faa:
        return len(set(re.findall(r'>[^_]+',faa.read())))

has_orthos = set([f[:-4] for f in os.listdir('correlatedEvolution/orthogroupFaa')])
puf_orthos = [n_orthologs(p.pp_id) for p in pufs if p.pp_id in has_orthos]
pkf_orthos = [n_orthologs(p.pp_id) for p in pkfs if p.pp_id in has_orthos]
fade = 0.5
fig,ax = plt.subplots()
ax.scatter([scatter(1) for _ in range(len(pkf_orthos))],
           pkf_orthos, 
           s = 1, color = 'k', marker = '.', alpha = fade)
ax.scatter([1 + scatter(1) for _ in range(len(puf_orthos))],
           puf_orthos, 
           s = 1, color = 'k', marker = '.', alpha = fade)
# ax.set_yscale('log')
ax.set_ylabel('Number of Orthologs')
ax.set_box_aspect(1)
ax.set_xticks((0,1),('PKFs','PUFs'))
fig.savefig('presentations_and_reports/paper_figures/nOrthologs.svg')

fig,ax = plt.subplots()
ax.scatter(np.quantile(pkf_orthos,np.linspace(.01,.99,100)),
           np.quantile(puf_orthos,np.linspace(.01,.99,100)),
           s = 1, color = 'w', marker = '.', zorder=2)
ax.scatter(np.quantile(pkf_orthos,.5),
           np.quantile(puf_orthos,.5),
           s = 3, color = 'r', marker = '.', zorder=2)
ax.text(np.quantile(pkf_orthos,.5) + 3, 
        np.quantile(puf_orthos,.5) -20, 
        'Median', fontsize = 6)
ax.plot((0,600),(0,600),'--r',linewidth = 0.5)
for i,line in enumerate(bootstrap_qq(pkf_orthos, puf_orthos)):
    ax.plot(line[0],line[1],'-k', linewidth = 0.2, alpha = 0.2, zorder=1)
ax.set_box_aspect(1)
ax.set_xlabel('PKF Orthologs')
ax.set_ylabel('PUF Orthologs')
ax.set_facecolor('lightgrey')
ax.set_title('Q-Q Plot of Ortholog Counts')
fig.savefig('presentations_and_reports/paper_figures/nOrthologs-QQ.svg')

######Figure S2A
sim_results = pd.read_csv('termCentric/SLT01_proteinSimilarityModelResults.tsv',sep = '\t')
sim_results = sim_results[sim_results['partition'] == 'Te']
true_scores = sim_results[sim_results['depth_long'] > 6]['predicted_depth']
false_scores = sim_results[sim_results['depth_long'] <= 6]['predicted_depth']


true_scores = SortedList(true_scores)
false_scores = SortedList(false_scores)
all_scores = true_scores + false_scores
tprs = [0]
fprs = [0]
n_true = len(true_scores)
n_false = len(false_scores)
for score in reversed(all_scores):
    fn = true_scores.bisect_left(score)
    tp = n_true - fn
    tn = false_scores.bisect_left(score)
    fp = n_false - tn
    
    tprs.append(TPR([tp,fp,tn,fn]))
    fprs.append(FPR([tp,fp,tn,fn]))
tprs.append(1)
fprs.append(1)

fig, ax = plt.subplots()
ax.plot(fprs,tprs,color = 'k', linewidth = .5, clip_on=False)
ax.set_box_aspect(1)
ax.plot((0,1),(0,1),'--r', linewidth = .5)
ax.set_ylabel('True Positive Rate')
ax.set_xlabel('False Positive Rate')
ax.set_ylim((0,1.005))
ax.set_xlim((-0.005,1))
ax.set_title('Functional Similarity Test Set ROC')
fig.savefig('presentations_and_reports/paper_figures/func_sim_ROC.svg')
print(f'AUC: {trapz(tprs,fprs)}')
