# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:59:59 2022

@author: Administrator
"""

import dill
import os
from collections import defaultdict
import re
import networkx as nx
import pandas as pd
# import math
from copy import deepcopy

#protein object holds data from proteins.dat, sourced from biocyc
class protein:
    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            return attr, value
    def __init__(self):
        self.minimal_annotation = []
        self.catalyzes = []
        self.features = []
        self.locations = ''
        self.comment = ''
        self.mol_weight_kd = 0
        self.mol_weight_seq = 0
        self.uni_id = ''
        self.predicted_annotation = []
        self.unannotated = True
        self.bc_id = ''
        self.name_ref = ''
    def parse_fasta(self, entry):
        self.pp_id = re.search(r'PP_\d{4}', entry).group()
        self.seq = entry.split('\n')[1].strip()
    #the biocyc file protein.dat is split on r'//' and this function takes a single resultant string to parse which represnets a protein entry
    def parse_biocyc(self, dat_entry):
        gostring = 'GO:\d{7}'
        fields = dat_entry.split('\n')
        for field in fields:
            if re.match('TYPES',field):
                self.types = re.findall('- .+',field)[0][2:]
            elif re.match('COMMENT',field):
                self.comment = re.findall('- .+',field)[0][2:]
            elif re.match('GO-TERMS',field):
                self.minimal_annotation.append(re.findall(gostring,field)[0])
                self.unannotated = False
            elif re.match('CATALYZES',field):
                self.catalyzes.append(re.findall('- .+',field)[0][2:])
            elif re.match('FEATURES',field):
                self.features.append(re.findall('- .+',field)[0][2:])
            elif re.match('LOCATIONS',field):
                self.locations = re.findall('- .+',field)[0][2:]
            elif re.match('MOLECULAR-WEIGHT-KD',field):
                self.mol_weight_kd = float(re.findall('- .+',field)[0][2:])
            elif re.match('MOLECULAR-WEIGHT-SEQ',field):
                self.mol_weight_seq = float(re.findall('- .+',field)[0][2:])
    def parse_uniprot(self, uni_entry):
        self.uni_id = uni_entry['Entry']
        self.bc_id = uni_entry['biocyc']
        self.pp_id = uni_entry['Gene names  (ordered locus )']
        self.uniprot_name = uni_entry['Entry name']
        self.name_ref = uni_entry['ref']
        self.name_gi = uni_entry['gi']
        self.name_primary = uni_entry['Gene names  (primary )']
        self.name_orf = uni_entry['Gene names  (ORF )']
        self.name_synonym = uni_entry['Gene names  (ORF )']
        self.name_protein = uni_entry['Protein names']
        self.seq = uni_entry['seq']
        if type(uni_entry['Gene ontology IDs']) == str:
            self.minimal_annotation.extend([t.strip() for t in uni_entry['Gene ontology IDs'].split(';')])
            self.unannotated = False
    #starting from the annotations given by biocyc this fills out the annotation list with all parent annotations and ensures that alt terms
    #are replaced with the primary synonym
    def full_annotation(self,GO,alt_dict):
        self.annotation = []
        for i,a in enumerate(self.minimal_annotation):
            if not a in list(GO.nodes()):
                if a in list(alt_dict.keys()):
                    self.minimal_annotation[i] = alt_dict[a]
                    self.annotation.append(nx.descendants(GO,alt_dict[a]))
                else:
                    self.minimal_annotation.remove(a)
            else:
                self.annotation.append(nx.descendants(GO,a))
        self.annotation = set([item for sublist in self.annotation for item in sublist])
        if any([go_dict[t].depth_long > 0 for t in self.annotation]):
            self.unannotated = False
        else:
            self.annotation = []
            self.unannotated = True
    def dump(self, file):
        dump_list = []        
        for attribute in self.__dict__.keys():
            dump_list.append('\t'.join([str(type(self.__dict__[attribute])), attribute, str(self.__dict__[attribute])]))
            with open(file,'w') as outfile:
                outfile.write('\n'.join(dump_list))

class term:
    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            return attr, value
    def __init__(self):
        self.alt_id = []
        self.subset = []
        self.synonym = []
        self.xref = []
        self.is_a = []
        self.intersection_of = []
        self.relationship = []
        self.consider = []
        self.id = ''
        self.name = ''
        self.namespace = ''
        self.definition = ''
        self.comment = ''
        self.info_content = 20
    #parses go.obo from the GO consortium, file must be split on '[Term]' and a single resultant string is passed as obo_entry
    def parse_input(self, obo_entry):
        gostring = 'GO:\d{7}'
        fields = obo_entry.split('\n')
        for field in fields:
            if re.match('id:',field):
                self.id = re.findall(gostring,field)[0]
            elif re.match('name:',field):
                self.name = re.findall(' .+',field)[0][1:]
            elif re.match('namespace:',field):
                self.namespace = re.findall(' .+',field)[0][1:]
            elif re.match('def:',field):
                self.definition = re.findall(' .+',field)[0][1:]
            elif re.match('comment:',field):
                self.comment = re.findall(' .+',field)[0][1:]
            elif re.match('alt_id:',field):
                self.alt_id.append(re.findall(gostring,field)[0])
            elif re.match('subset:',field):
                self.subset.append(re.findall(' .+',field)[0][1:])
            elif re.match('synonym:',field):
                self.synonym.append(re.findall(' .+',field)[0][1:])
            elif re.match('xref:',field):
                self.xref.append(re.findall(' .+',field)[0][1:])
            elif re.match('is_a:',field):
                self.is_a.append(re.findall(gostring,field)[0])
            elif re.match('intersection_of:',field):
                self.intersection_of.append(re.findall(' .+',field)[0][1:])
            elif re.match('relationship:',field):
                self.relationship.append(re.findall(' .+',field)[0][1:])
            elif re.match('consider:',field):
                self.consider.append(re.findall(gostring,field)[0])


os.chdir('C:/Users/Administrator/Documents/Hettich/SLT01')

############################## GO graph setup ###########################
#parses go.obo into list of go term objects and makes a dictionary of term ids : object
with open('baseline_information/go.obo', 'r') as obo:
    termlist = obo.read().split('[Term]')[1:]
del obo
terms = []
termlist[-1] = termlist[-1].split('[Typedef]')[0]
for i,t in enumerate(termlist):
    if not re.search('is_obsolete: true', t):
        terms.append(term())
        terms[-1].parse_input(termlist[i])
del termlist

#sets up the GO graph as a networkx DiGraph with GO term IDs for nodes. only is_a relationships are considered
term_edges = []
for t in terms:
    for a in t.is_a:
        term_edges.append(tuple([t.id,a,{type:'is_a'}]))
go_g = nx.DiGraph(term_edges)
del term_edges

#translates between alternative IDs for go terms and their primary synonym
alt_id_dict = {}
for t in terms:
    for alt in t.alt_id:
        alt_id_dict[alt] = t.id

#this adds a depth value for each term in the GO which is defined as the shortest path length to its namespace
top_terms = [t for t in terms if go_g.out_degree(t.id) == 0]
ontology = {'biological_process':'GO:0008150',
            'cellular_component':'GO:0005575',
            'molecular_function':'GO:0003674'}
for t in terms:
    if t in top_terms:
        t.depth_short = 0
        t.depth_long = 0
    else:
        t.depth_short = nx.shortest_path_length(go_g, t.id, ontology[t.namespace])
        t.depth_long = max([len(p) for p in nx.all_simple_paths(go_g, t.id, ontology[t.namespace])]) - 1

go_dict = {t.id:t for t in terms}

# # calculates the infromation content of GO terms from how frequently they are annotated to bacterial proteins in swissprot
# alts = set(alt_id_dict.keys())
# go_freqs = defaultdict(lambda: 0)
# prot_counter = 0
# with open('other_files/uniprot-swissProt_GOannotations.tab','r') as uniprot_file:
#     for line in uniprot_file:
#         if 'Bacteria' in line:
#             min_annotation = [alt_id_dict[a] if a in alts else a for a in re.findall(r'GO:\d{7}', line)]
#             annotation = set([t for term in min_annotation if term in go_g.nodes for t in nx.descendants(go_g,term)])
#             for term in annotation:
#                 go_freqs[term] += 1
#             if annotation:
#                 prot_counter += 1
# with open('other_files/uniprot-trEMBL_GOannotations.tab','r') as uniprot_file:
#     for line in uniprot_file:
#         if 'Bacteria' in line:
#             min_annotation = [alt_id_dict[a] if a in alts else a for a in re.findall(r'GO:\d{7}', line)]
#             annotation = set([t for term in min_annotation if term in go_g.nodes for t in nx.descendants(go_g,term)])
#             for term in annotation:
#                 go_freqs[term] += 1
#             if annotation:
#                 prot_counter += 1

# for term in terms:
#     term.info_content = -math.log2(go_freqs[term.id]/prot_counter) if go_freqs[term.id] > 0 else -math.log2(1/prot_counter)

# max_info = max([t.info_content for t in terms])
# max_depth = max([t.depth_long for t in terms])
# for term in terms:
#     term.lhalf_sim = (math.sqrt(term.info_content/max_info) + math.sqrt(term.depth_long/max_depth))**2

# max_lhalf = max([t.lhalf_sim for t in terms])
# for term in terms:
#     term.lhalf_sim = int(round((term.lhalf_sim/max_lhalf)*max_depth))

# del term, go_freqs, prot_counter, a, alt, i, t, line, uniprot_file, alts, min_annotation, annotation

############################## protein object list setup ###########################
#sets up the list of protein objects that holds all protein data
proteins = []
with open('baseline_information/Pseudomonas_putida_KT2440_110.faa','r') as faa:
    entries = faa.read().split('>')[1:]
for entry in entries:
    proteins.append(protein())
    proteins[-1].parse_fasta(entry)
fasta_pps = [p.pp_id for p in proteins]

with open('baseline_information/Pseudomonas_putida_KT2440_7360.faa','r') as faa:
    entries = faa.read().split('>')[1:]
for entry in entries:
    if re.search(r'PP_\d{4}',entry).group() not in fasta_pps:
        proteins.append(protein())
        proteins[-1].parse_fasta(entry)
fasta_pps = [p.pp_id for p in proteins]

names_in = pd.read_csv('baseline_information/names_seq_input.csv', dtype=str)
for i in range(len(names_in['Entry'])):
    entry = names_in.iloc[i,:]
    #uniprot collapses some duplicated genes into a single entry
    for pp in entry[3].split(' '):
        e = deepcopy(entry)
        e[3] = pp
        if pp in fasta_pps:
            next(p for p in proteins if p.pp_id == pp).parse_uniprot(e)
        else:
            proteins.append(protein())
            proteins[-1].parse_uniprot(e)
del fasta_pps, faa

#dicts for translating from other names to uniprot IDs, uniprot IDs will be the internal naming scheme for the whole script
list_of_uni = names_in['Entry'].to_list()
pp_uni = {}
for i,loci in enumerate(names_in['Gene names  (ordered locus )'].to_list()):
    if loci.split() == 1:
        pp_uni[loci] = list_of_uni[i]
    else:
        for locus in loci.split():
            pp_uni[locus] = list_of_uni[i]
del names_in, list_of_uni
#biocyc IDs to uniprot IDs
bc_uni = {bc:uni for bc,uni in zip([b.bc_id for b in proteins],[u.uni_id for u in proteins])}
#refseq IDs to uniprot IDs
ref_uni = {ref:uni for ref,uni in zip([r.name_ref for r in proteins],[u.uni_id for u in proteins])}
#associates uniprot IDs with their protein object to speed up searches
uni_dict = {uni:[pr for pr in proteins if pr.uni_id == uni][0] for uni in [u.uni_id for u in proteins]}

#parse protein.dat into a list of protein objects
#some entries do not have a uniprot ID so they must be matched to a protein object through 
with open('baseline_information/proteins.dat','r') as dat:
    protlist = dat.read().split(r'//')
del dat
gene_pattern = r"GENE - .+"
complex_pattern = r'UNIQUE-ID - CPLX'
for p in protlist:
    if re.search(gene_pattern, p):
        gene_name = re.search(gene_pattern, p).group()[7:]
    else:
        gene_name = ''
    if re.search(r'UNIPROT ".+"', p):
        uni_name = re.search(r'UNIPROT ".+"', p).group()[9:-1]
    else:
        uni_name = ''
    if re.search(r'"NP_.+"',p):
        ref_name = re.search(r'"NP_.+"',p).group()[1:-1]
    else:
        ref_name = ''
    if gene_name in [pr.bc_id for pr in proteins]:
        prot = uni_dict[bc_uni[gene_name]]
    elif uni_name in [pr.uni_id for pr in proteins]:
        prot = uni_dict[uni_name]
    elif ref_name in [pr.name_ref for pr in proteins]:
        prot = uni_dict[ref_uni[ref_name]]
    elif re.search(complex_pattern, p) is not None:
        component_list = re.findall(r'COMPONENTS - .{5}-\d+-', p)
        go_list = re.findall(r'GO:\d{7}', p)
        if type(go_list) == str:
            go_list = [go_list]
        for component in component_list:
            component = component[13:-1]
            prot = uni_dict[bc_uni[component]]
            prot.minimal_annotation.extend(go_list)
        continue
    else:
        print(p)
        continue
    prot.parse_biocyc(p)
# del protlist, uni_name, ref_name, gene_name, gene_pattern, complex_pattern, go_list

#GO annotations from Pseudomonas Genome DB. The union of this, proteins.dat from BioCyc, and annotations in uniprot will be the initial genome annotation
pgdb_go = pd.read_csv('baseline_information/gene_ontology_csv.csv')
pgdb_go['uni'] = [pp_uni[p] for p in pgdb_go['Locus Tag']]
uni_names = [pr.uni_id for pr in proteins]
for i in range(len(pgdb_go['Locus Tag'])):
    if  pgdb_go.iloc[i,12] in uni_names:
        prot = uni_dict[pgdb_go.iloc[i,12]]
        prot.minimal_annotation.append(pgdb_go.iloc[i,4])
    else:
        print(pgdb_go.iloc[i,12])
del uni_names, pgdb_go

#GO annotations from NetGo
uni_prot = defaultdict(lambda: [])
for p in proteins:
    uni_prot[p.uni_id].append(p)
for result_file in [f for f in os.listdir('baseline_information') if f.startswith('result_')]:
    with open('baseline_information/'+result_file) as netgo:
        for line in netgo:
            entries = line.split('\t')
            if len(entries) > 1:
                if float(entries[2]) > .90:
                    for p in uni_prot[entries[0]]:
                        p.minimal_annotation.append(entries[1])

#full_annotations makes sure that all GO accessions are translated to the main synonym, that there are no repeated terms,
#and that all parent terms of an included term are also included
for prot in proteins:
    prot.full_annotation(go_g, alt_id_dict)

# del p, i, locus, loci, entry, prot, component_list, pp, e

with open('pickled_objects/proteins.dill','wb') as picklefile:
    dill.dump(proteins,picklefile)
with open('pickled_objects/terms.dill','wb') as picklefile:
    dill.dump(terms,picklefile)
with open('pickled_objects/go_g.dill','wb') as picklefile:
    dill.dump(go_g,picklefile)


