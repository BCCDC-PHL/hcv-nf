#!/usr/bin/env python3

import sys as sys
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd
import argparse

def main(args):
    
    bedfile = pd.read_csv(args.bed,sep="\t",header = None)
   

    with open(args.name) as f:
        genotype_name = f.readlines()
        
    genotype_name = [item.strip() for item in genotype_name]
    segment_name = [item.split('|')[2] + '_'+ item.strip().split('|')[1] for item in genotype_name]
    sample_name = genotype_name[0].split('_')[0]

    bed = bedfile[bedfile[0].isin(segment_name)]

    res = {segment_name[i]: genotype_name[i] for i in range(len(segment_name))}
    bed.replace({0: res},inplace=True)

    bed.to_csv(sample_name+'_amplicon.bed', sep='\t',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed', required=True)
    parser.add_argument('--name', required=True)

    args = parser.parse_args()
    main(args)
