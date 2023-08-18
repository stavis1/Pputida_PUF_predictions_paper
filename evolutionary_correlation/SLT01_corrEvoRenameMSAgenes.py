import sys
import os
import re
from multiprocessing import Pool
from multiprocessing import cpu_count

def rename_genes(msa_filename):
    with open(msa_dir + msa_filename,'r') as msa_in:
        with open(out_dir + msa_filename, 'w') as msa_out:
            for line in msa_in:
                if line.startswith('>'):
                    gene = gene = re.match(r'>\w+', line).group()[1:]
                    line = '>' + gene_species_dict[gene] + '\n'
                msa_out.write(line)
    return

msa_dir = sys.argv[1]
faa_dir = sys.argv[2]
out_dir = sys.argv[3]
if len(sys.argv) > 4:
    cores = int(sys.argv[4])
else:
    cores = cpu_count()

msa_files = [f for f in os.listdir(msa_dir) if f[-3:] == 'faa']
faa_files = [f for f in os.listdir(faa_dir) if f[-3:] == 'faa']

gene_species_dict = {}
for faa_file in faa_files:
    with open(faa_dir + faa_file, 'r') as faa:
        for line in faa:
            if line.startswith('>'):
                gene = re.match(r'>\w+', line).group()[1:]
                gene_species_dict[gene] = faa_file[:-4]

with Pool(cores) as p:
    p.map(rename_genes, msa_files)


