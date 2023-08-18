# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:33:23 2022

@author: Administrator
"""
import os
import sys
import pandas as pd
import requests
import re
from multiprocessing import Pool
import time
from requests.adapters import HTTPAdapter, Retry

def uniprot_map(i):
    print(f'starting uniprot id lookup for set {i}')
    ids = id_sets[i]
    
    #initiate the ID mapping job on uniprot's servers
    map_job = requests.post("https://rest.uniprot.org/idmapping/run",
                            data={"from":"PDB", "to":"UniProtKB", "ids":",".join(ids)})
    job_id = map_job.json()['jobId']
    
    #wait until job completes
    while True:
        request = requests.get(f"https://rest.uniprot.org/idmapping/status/{job_id}")
        j = request.json()
        if "jobStatus" in j.keys():
            if j["jobStatus"] == 'RUNNING':
                time.sleep(interval)
            else:
                break
        else:
            break
    
    #retrieve results and retry if there's an issue
    results = requests.get(f"https://rest.uniprot.org/idmapping/results/{job_id}?size={size}").json()
    if not 'results' in results.keys():
        print('BAD BATCH:')
        print(results)
        results = requests.get(f"https://rest.uniprot.org/idmapping/results/{job_id}?size={size}").json()
    
    #parse output
    id_map = {}
    for result in results['results']:
        pdb = result['from']
        uniprot = result['to']
        id_map[pdb] = uniprot
    mapped_ids = set(id_map.keys())
    unmapped = [i for i in ids if i not in mapped_ids]
    print('number of ids mapped: {n}'.format(n=len(mapped_ids)))
    return (id_map, unmapped)
    

def pdb_lookup(chain):
    #construct URL
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=5)
    session.mount('https://', adapter)
    SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"
    UNIPROT = "/mappings/uniprot"
    full_url = "%s/%s/%s" % (SERVER_URL, UNIPROT, chain[:4])
    
    #attempt to retrieve data, if there's no uniprot ID in the output or if the model ID is bad return None
    try:
        json_results = session.get(full_url, timeout=10).json()
        results = json_results[list(json_results.keys())[0]]['UniProt']
    except Exception as e:
        print('bad chain : {c}\n{e}'.format(c=chain,e=e))
        return None
    
    #look for the specific chain in the results, if it can't be found return None
    unis = list(results.keys())
    for u in unis:
        if re.search('\'chain_id\': \'{c}'.format(c=chain[5:].lower()),str(results[u]).lower()):
            return u
    return None

def get_data(protein):
    #retrieve all uniprot data on a particular protein and extract the sequene, lineage, and GO terms
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=5)
    session.mount('https://', adapter)
    data = session.get(f'https://www.uniprot.org/uniprotkb/{protein}.xml?include=yes').text
    seq = re.search(r'([^>]+)</sequence>',data).group(1)
    go_terms = ';'.join(re.findall(r'GO:\d{7}',data))
    lineage = ';'.join(re.findall('<taxon>([^<]+)<',data))
    return (seq,go_terms,lineage)



#######
interval = 10
size = 500
max_uni = 100000
if len(sys.argv) > 4:
    cores = int(sys.argv[4])
else:
    cores = os.cpu_count()
if len(sys.argv) > 5:
    pdb_cores = int(sys.argv[5])
else:
    pdb_cores = cores
retries = Retry(total=5,
                backoff_factor=0.1,
                status_forcelist=[ 500, 502, 503, 504 ])
#######

puf_dir = os.path.abspath(sys.argv[1]) + '/'
ann_dir = os.path.abspath(sys.argv[2]) + '/'
out_dir = os.path.abspath(sys.argv[3]) + '/'

puf_files = [puf_dir + f for f in os.listdir(puf_dir) if f.endswith('.txt')]
ann_files = [ann_dir + f for f in os.listdir(ann_dir) if f.endswith('.txt')]

#extract all chain IDs from RUPEE hit files
ids = []
for puf in puf_files:
    hits = pd.read_csv(puf)
    ids.extend([c[:4].upper() + ':' + c[4:].upper() for c in hits['Chain Id']])
for ann in ann_files:
    hits = pd.read_csv(ann)
    ids.extend([c[:4].upper() + ':' + c[4:].upper() for c in hits['Chain Id']])
ids = set(ids)

#the uniprot streaming node can't deal with batches larger than 500
#so break up the searches into chunks
n_sets = 1 + int(len(ids)/size)
id_sets = [[p for j,p in enumerate(ids) if j%n_sets == i] for i in range(n_sets)]
print('there are {n} id sets'.format(n=n_sets))

#search for ID mappings using uniprot
with Pool(cores) as p:
    uni_out = p.map(uniprot_map,range(n_sets))

#pool the results
map_accumulator = {}
unmapped_accumulator = []
for result in uni_out:
    map_accumulator.update(result[0])
    unmapped_accumulator.extend(result[1])

print('total number of mapped ids: {n}'.format(n=len(map_accumulator.keys())))
print('total number of unmapped ids sent to pdb lookup: {n}'.format(n=len(unmapped_accumulator)))

#search for ID mappings on PDB directly (this is slower so it's done only on the unmapped IDs)
with Pool(pdb_cores) as p:
    pdb_results = p.map(pdb_lookup, unmapped_accumulator)

#pool results
map_accumulator.update({pdb:uni for pdb,uni in zip(unmapped_accumulator,pdb_results) if uni is not None})

#chain_map.tsv only holds chain ID -> uniprot ID
with open(out_dir + 'chain_map.tsv','w') as map_file:
    map_file.write('pdb\tuniprot\n')
    map_file.write('\n'.join([k + '\t' + map_accumulator[k] for k in map_accumulator.keys()]))

#the pairwise jobs file is needed for running pairwise sequence alignments
#on all mappable hits
pairwise_jobs = []
mapped_chains = set(map_accumulator.keys())
#collect jobs coming from PUF hits
for puf in puf_files:
    prot = re.search(r'PP_\d{4}',puf).group()
    hits = pd.read_csv(puf)
    hits['Chain Id'] = [c[:4].upper() + ':' + c[4:].upper() for c in hits['Chain Id']]
    hits = hits[[c in mapped_chains for c in hits['Chain Id']]]
    pairwise_jobs.extend(['{p}\t{u}'.format(p=prot,u=map_accumulator[chain]) for chain in hits['Chain Id']])
#collect jobs coming from annotated hits
for ann in ann_files:
    prot = re.search(r'PP_\d{4}',ann).group()
    hits = pd.read_csv(ann)
    hits['Chain Id'] = [c[:4].upper() + ':' + c[4:].upper() for c in hits['Chain Id']]
    hits = hits[[c in mapped_chains for c in hits['Chain Id']]]
    pairwise_jobs.extend(['{p}\t{u}'.format(p=prot,u=map_accumulator[chain]) for chain in hits['Chain Id']])
#save all jobs to a single file for SLT01_RUPEEpairwiseAlignments.py
with open(out_dir + 'pairwise_jobs.txt','w') as jobs_file:
    jobs_file.write('\n'.join(set(pairwise_jobs)))

#get data on all mapped uniprot hits
unis = list(set(map_accumulator.values()))
with Pool(cores) as p:
    uniprot_data = p.map(get_data,unis)

#fasta file for SLT01_RUPEEpairwiseAlignments.py
with open(f'{out_dir}uniprot_hits.fasta','w') as fasta_file:
    fasta_file.write('\n'.join(['>{u}\n{s}'.format(u=u,s=d[0]) for u,d in zip(unis,uniprot_data)]))

#lineage and GO data for SLT01_RUPEEsemisupervisedInput.py
with open(f'{out_dir}uniprot_hits.tsv','w') as data_file:
    data_file.write('uniprot\tGO\tlineage\n')
    data_file.write('\n'.join(['{u}\t{g}\t{l}'.format(u=u,g=d[1],l=d[2]) for u,d in zip(unis,uniprot_data)]))


