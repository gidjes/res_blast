import argparse
import os.path
import glob
import shutil
import subprocess 
from pathlib import Path
import pandas as pd
import csv

# Idea is to extract the contigs where blast result is found for in the .csv blast results Optional run mash straight after
blastn_header = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore', "ab_class"]
blastdir = os.path.abspath("./output/tmp")
outdir = os.path.abspath('./output')
all_files = glob.glob(f"{blastdir}/*/*.csv")

def csv_to_list(file):
    with open(file) as f:
        if os.path.getsize(file) > 0:
            reader = csv.reader(f)
            ab_class = os.path.basename(os.path.dirname(file))
            full_data = [row+[ab_class] for row in reader]
            return full_data

def csv_to_list_single(file):
    with open(file) as f:
        if os.path.getsize(file) > 0:
            reader = csv.reader(f)
            data = list(reader)[0]
            ab_class = os.path.basename(os.path.dirname(file))
            data += [ab_class]
            return data


data_single = [csv_to_list_single(f) for f in all_files if csv_to_list_single(f) != None]
data_full_all = [csv_to_list(f) for f in all_files if csv_to_list(f) != None]
data_full = [item for sublist in data_full_all for item in sublist]

df = pd.DataFrame(data_full, columns=blastn_header)
print(df)
df_single = pd.DataFrame(data_single, columns=blastn_header)

contig_file = f"{outdir}/contig_found.txt"
uni_contig_file = f"{outdir}/unique_contig_found.txt"
all_100_hits = f"{outdir}/all_100_contig_found.csv"
df.pident = df.pident.astype(float)
df_single.pident = df_single.pident.astype(float)
df = df[df.pident > 100.000]
df_single = df_single[df_single.pident > 100.000]

uni_list = list(set(df['qseqid'].tolist()))
uni_tn_list = list(set([a.split('_')[0] for a in uni_list]))
print(uni_tn_list)
df_uni = pd.DataFrame(uni_list)

# saving the dataframe s
df_uni.to_csv(uni_contig_file, index=None) 
df['qseqid'].to_csv(contig_file, index=None)
df.to_csv(all_100_hits, index=None)
typened_list = [e.split('_')[0] for e in df['qseqid'].tolist()]
all_hybrid_files = glob.glob(f"./hybrid_fasta/*")

destination_paths = [f"{outdir}/{b}_extract.fasta" for b in set(typened_list)]

def create_hybrid_list(inputfiles, listfound):
    temp_all_hybrid_files_found = []
    temp_typened_file_dict = {}
    for f in inputfiles:
        searchstr = os.path.basename(f).split('_')
        for x in searchstr:
            if x in listfound:
                temp_all_hybrid_files_found.append(f)
                temp_typened_file_dict[f] = x

    return temp_all_hybrid_files_found, temp_typened_file_dict

create_hybrid_list_output = create_hybrid_list(all_hybrid_files, typened_list)
hybrid_found_list = create_hybrid_list_output[0] # This list contains all 'key_contig' entries to extract to a new file using seqtk subseq.
typened_file_dict = create_hybrid_list_output[1] # Will have the correct type-ned key associated with the file name used

# print(contig_file)

def gen_meta_contig_file(flagcontigkey):
    fasta_contig = os.path.abspath(f"{outdir}/fasta_contig_found")
    Path(os.path.abspath(fasta_contig)).mkdir(parents=True, exist_ok=True)
    for idx, file in enumerate(hybrid_found_list):
        savename = f"{fasta_contig}/{typened_file_dict[file]}_blasthits.fasta"
        # if savename.exists():
        subprocess_cmd = f"seqtk subseq {file} {flagcontigkey} >> {savename}"
        subprocess.Popen(subprocess_cmd, shell=True, stdout=subprocess.PIPE)

gen_meta_contig_file(uni_contig_file)