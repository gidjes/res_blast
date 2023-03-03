import argparse
import os.path
import glob
import shutil
import subprocess 
import pandas as pd
# import pymssql
# import pathlib
from pathlib import Path

# Example output:
# query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 22006156_2,blaNDM-1_1_FN396876,100.000,813,0,0,89183,89995,1,813,0.0,1502
# Example usage
# python scripts/blastn2sequence.py

out = os.path.abspath("./output/tmp/")
fastas = Path(os.path.abspath("../fastas"))
list_of_files = glob.glob(os.path.abspath(f"{fastas}/*"))

resfinder_path = os.path.abspath("./resfinder_db")
search_string = resfinder_path + "/*.fsa"
files = glob.glob(search_string)

for x in files:
    ab_class = x.split("/")[-1].replace(".fsa", "")
    output_path = out + "/" + ab_class
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    for key in list_of_files:
        import subprocess
        basename = key.split('/')[-1].split('.')[0]
        outputname = f"{output_path}/{basename}_blastn.csv"
        subprocess_cmd_blast2input = f"bsub -q bio ./logs/Process_log.txt -e ./logs/Error_log.txt\
                -n 1 -R 'rusage[mem=12G]' -R 'span[hosts=1]' -W 15 \
                \"blastn -query {key} \
                -subject {x} \
                -out {outputname} \
                -outfmt 10\""
        subprocess.Popen(subprocess_cmd_blast2input, shell=True, stdout=subprocess.PIPE)
    print(f"Jobs sent to cluster for {len(list_of_files)} isolates and {ab_class}.")
