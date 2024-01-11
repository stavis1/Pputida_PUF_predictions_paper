# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:58:42 2022

@author: Administrator
"""
from collections import defaultdict
import pandas as pd
import itertools
import math
import os
import re
import numpy as np
from numba import jit
import dill

os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01')

#################################################
#expression correlation given by the mean of the set of pearson correlation coefficents of protein expression z scores across three proteomics datasets
#and jaccard similarity of binary presence/absence vector across conditions

experiments = ['PXD013011', 'PXD016028', 'PXD016114']
proteins = set([])
for experiment in experiments:
    results = pd.read_csv('SLT01_coexpression/SLT01_flashLFQ/' + experiment + '_QuantifiedProteins.tsv', sep = '\t', index_col = 'Protein Groups')
    proteins.update(results.index)
#this protein list will be used to give a consistant ordering and index to protein data, which is necessary as not all proteins are observed in all experiments
proteins = list(proteins)

corr_data = []
presence_data = pd.DataFrame(index = proteins)
for experiment in experiments:
    #QuantifiedProteins.tsv is the raw output of flashLFQ with no bayesian fold change analysis
    #see the .toml files in the same directory for details
    results = pd.read_csv('SLT01_coexpression/SLT01_flashLFQ/' + experiment + '_QuantifiedProteins.tsv', sep = '\t', index_col = 'Protein Groups')
    results = results.drop([c for c in results.columns if not c.startswith('Intensity')],axis = 1)
    conds = set([re.search(r'C(\d{1,2})R', c).group(1)for c in results.columns])
    combined = pd.DataFrame(index = proteins)
    #intensity is averaged for each condition
    #proteins not observed in the experiment are given zero intensities for all conditions this is necessary for combining experiment data
    minval = min([v for v in results.to_numpy().flatten() if v > 0.0])
    for c in results.columns:
        results[c] = np.log([v if v > 0 else minval for v in results[c]])
    for c in conds:
        cond_frame = results.filter([col for col in results.columns if re.search(r'C' + c + 'R', col)])
        combined['C' + c] = pd.concat([cond_frame.mean(axis = 1),pd.Series({p:minval for p in proteins if p not in results.index})])
    #presenece holds binary vectors indicating whether a protein was observed in a condition
    minval = min([v for v in results.to_numpy().flatten() if v > 0.0])
    presence = combined.loc[proteins].apply(lambda x: [bool(i > minval) for i in x], result_type= 'broadcast').astype(bool)
    presence.columns = [c + experiment for c in presence.columns]
    presence_data = pd.concat([presence_data,presence],axis = 1)
    #transposes the combined dataframe then gives all pairwise correlations between columns 
    #corr_matrix is a square numpy matrix with column and row indices correspoding to the protiens list
    corr_matrix = combined.loc[proteins].T.corr(method = 'spearman').to_numpy()
    #unobserved proteins (i.e. all 0) have no variance so they return nan, these are replaced with zero correlation
    corr_matrix = np.where(np.isnan(corr_matrix),0,corr_matrix)
    corr_data.append(corr_matrix)

#takes the jaccard similarity of the presence/absence vectors for each protien pair
#@jit is a decorator to allow just in time compilation by numba
#if this doesn't work you can comment the @jit line out and the code should still run, just slower
@jit(nopython = True, parallel = True, fastmath = True)
def np_jaccard(x,y):
    x = np.array([bool(i) for i in x])
    y = np.array([bool(i) for i in y])
    if np.double(np.bitwise_or(x, y).sum()) == 0:
        return 0
    else:
        return np.double(np.bitwise_and(x, y).sum()) / np.double(np.bitwise_or(x, y).sum())

ex_jaccard = presence_data.T.corr(np_jaccard).to_numpy()

#the list of matricies is stacked into a three dimensional matrix with axis 2 representing experiment
#the data are combined using the square of the sum of the square roots of the correlation coefficent magnitudes. 
#The signs of each experiment are propagated to the square roots so that the sum of sqrts respects concordance. 
#The overall sign is given by the sign of the sum.
#this is a modified l1/2 norm which, unlike the l2 norm (euclidian distance), has the effect of suppressing discordant values.
#the modified part of this norm is how it handles negative numbers which would otherwise give difficult to interpret complex outputs or NaNs
# @jit(nopython = True, parallel = True, fastmath = True)
# def mod_lhalf(x):
#     return np.square(np.sum(np.sqrt(np.abs(x)) * np.sign(x))) * np.sign(np.sum(x))

corrs = np.apply_along_axis(np.mean, axis = 2, arr = np.stack(corr_data, axis = 2))

with open('edgelists/SLT01_coEx-edgelist.txt','w') as edgelist:
    edgelist.write('\t'.join(['protein1','protein2','coEx','exJaccard']) + '\n')
    for pair in itertools.combinations(range(len(proteins)), 2):
        edgelist.write('\t'.join([proteins[pair[0]],proteins[pair[1]],str(corrs[pair]),str(ex_jaccard[pair])]) + '\n')

#################################################
#evolutionary correlations measured as the size of the intersection of taxa with identified orthologs
#and the Robinson-Foulds weighted cluster metric on the gene trees pruned down to the intersection taxa set

#formats the output of SLT01_corrEvoTreeCMP-fullComparison.py into an edgelist
with open('correlatedEvolution/SLT01_corrEvo-rawRCW.tsv','r') as infile:
    infile.readline()
    with open('edgelists/SLT01_corrEvo-edgelist.txt','w') as edgelist:
        edgelist.write('\t'.join(['protein1','protein2','shared_taxa','rcw']) + '\n')
        for line in infile:
            elms = line.split('\t')
            edgelist.write('\t'.join([elms[1][:7],elms[2][:7],elms[5],elms[6].strip()]) + '\n')

#################################################
#protein sequence similarity calculated by a all-v-all pairwise diamond search of the P. putida proteome

#formats the raw diamond output
with open('SLT01_seqSimilar/SLT01_seqSimilar-diamond.txt','r') as dmnd:
    with open('edgelists/SLT01_seqSimilar-edgelist.txt', 'w') as edgelist:
        edgelist.write('\t'.join(['protein1','protein2','bitscore']) + '\n')
        for line in dmnd:
            if '*' in line:
                continue
            edge = []
            #protein names (locus codes)
            edge.extend(re.findall(r'PP_\d{4}', line))
            # #e-value and bit score
            edge.append(re.findall(r'\d+\.[\de+-]+',line)[-1])
            #calculates the shannon information content from the evalue so that the numbers are easier for the BART model to work with
            #some e-values are 0 due to float inprecision which would give infinity bits, these are instead given the max observed information content
            # edge[2] = str(-math.log2(float(edge[2]))) if not float(float(edge[2])) == 0 else '1044' #this was the max value observed
            edgelist.write('\t'.join(edge) + '\n')

#################################################
#operon co-membership from operon annotations given by rockhopper

def rename(gene):
    #rockhopper uses gene symbols where it can but locus codes are the lingua franca of my analysis
    if gene.startswith('PP_'):
        return gene
    else:
        return gene_dict[gene]

gene_dict = {}
#this file contains both locus codes and gene symbols
with open('SLT01_rockhopperOperonAnalysis/Rockhopper_Results/genomes/Pseudomonas_putida_KT2440/NC_002947.ptt','r') as ptt:
    for line in ptt:
        entries= line.split('\t')
        if len(entries) > 5:
            if entries[4] != '-':
                gene_dict[entries[4]] = entries[5]
        

#this is the raw rockhopper output
with open('SLT01_rockhopperOperonAnalysis/Rockhopper_Results/NC_002947_operons.txt','r') as op_file:
    with open('edgelists/SLT01_operon-edgelist.txt','w') as edgelist:
        edgelist.write('\t'.join(['protein1','protein2','operon']) + '\n')
        op_file.readline()
        for line in op_file:
            genes = line.split('\t')[-1].strip().split(', ')
            genes = [rename(gene) for gene in genes]
            # genes = [gene for gene in genes if gene is not None]
            for pair in itertools.combinations(genes, 2):
                edgelist.write('\t'.join(pair + tuple(['1'])) + '\n')
            
#################################################
#categorical similarity score using the jaccard distance of InterProScan features

def update_signalp_data(row):
    global cat_data
    #row[2] is either signal_peptide or lipoprotein_signal_peptid
    #row[8] indicates if signal_peptide is a TAT peptide
    cat_data[row[0]].add(row[2] + row[8])

def update_interpro_data(row):
    global cat_data
    cat_data[row['Locus Tag']].add(row['External Signature Accession'])

def jaccard_d(comp):
    #jaccard distance for a vector of weighted elements is 1-(sum(intersection)/sum(union))
    #here the element weights are the shannon information content of the feature
    return 1 - sum([features[f] for f in cat_data[comp[0]].intersection(cat_data[comp[1]])])/sum([features[f] for f in cat_data[comp[0]].union(cat_data[comp[1]])])

signalp = pd.read_csv('jaccard/Pseudomonas_putida_KT2440_110.gff3', sep = '\t', header=None, skiprows=1)
interpro = pd.read_csv('jaccard/features.csv')
cat_data = defaultdict(lambda: set([]))
signalp.apply(update_signalp_data, axis = 1)
interpro.apply(update_interpro_data, axis = 1)
#this produces a table of all observed sequence features and their frequencies
features = pd.Series([elm for loc in cat_data.keys() for elm in cat_data[loc]]).value_counts().to_dict()
n_prots = len(cat_data.keys())
#the shannon information content of each feature is calculated based on its observation frequency in the data
for f in features.keys():
    features[f] = -math.log2(features[f]/n_prots)
with open('edgelists/SLT01_jaccardDistanceLayer-edgelist.txt','w') as edgelist:
    edgelist.write('\t'.join(['protein1','protein2','jaccard']) + '\n')
    for comp in itertools.combinations(cat_data.keys(), 2):
        edgelist.write('\t'.join([comp[0],comp[1],str(jaccard_d(comp))]) + '\n')

#################################################
#STRINGdb data

str_db = pd.read_csv('STRINGdb_data/160488.protein.links.full.v11.5.txt', sep = ' ')
str_db['protein2'] = [p[7:] for p in str_db['protein2']]
str_db['protein1'] = [p[7:] for p in str_db['protein1']]
str_db[str_db.columns[2:]] = str_db[str_db.columns[2:]]/1000
str_db = str_db.drop(['coexpression','experiments'],axis= 1)
str_db.to_csv('edgelists/SLT01_stringDBfull-edgelist.txt', sep = '\t', index = False)

#################################################
#TM align output from alphafold predictions

# with open('pickled_objects/proteins.dill','rb') as pickle:
#     proteins = dill.load(pickle)

# uni_pp = {p.uni_id:p.pp_id for p in proteins}

tm_data = pd.read_csv('alphafold/trimmed_collected_out.txt', sep = '\t')
tm_data['#PDBchain1'] = [re.search('PP_\d{4}',f).group() for f in tm_data['#PDBchain1']]
tm_data['PDBchain2'] = [re.search('PP_\d{4}',f).group() for f in tm_data['PDBchain2']]
tm_data['max_tm'] = [max([t1,t2]) for t1,t2 in zip(tm_data['TM1'],tm_data['TM2'])]

edgelist = tm_data[['#PDBchain1','PDBchain2','max_tm','RMSD']]
edgelist.columns = ['protein1','protein2','max_tm','rmsd']
edgelist.to_csv('edgelists/SLT01_alphafold-edgelist.txt', sep = '\t', index = False)




