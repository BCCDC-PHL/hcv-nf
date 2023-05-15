#!/usr/bin/env python3

import argparse
import pandas as pd 
import numpy as np

def main(args):
    #cols = 'subject_taxids,lineage,name,rank'.split(',')
    #taxon_results = pd.read_csv('F1910235_taxon_results.txt',sep = '\t', names=cols)
    #taxon_results = pd.read_csv(args.taxonresult,sep = '\t', names=cols)

    type = {}
    with open(args.taxonresult, 'r') as f:
        for line in f.readlines():
            subject_id = line.split('\t')[0]
            lineages = line.split('\t')[1].split(';')
            if len(lineages) > 10 :
                lineage = lineages[9]
            else:
                lineage = lineages[-1].replace(',','')

            type[subject_id] = lineage

    seq_description = {}
    with open(args.descriptionfile, 'r') as f:
        for line in f.readlines():
            subject_accession = line.split(' ')[0].replace('>','')

            description = line.split(' ')[1::]
            description = ' '.join(description)
            seq_description[subject_accession] = description

    blast_results = pd.read_csv(args.blastresult)
    blast_results = blast_results.assign(subject_taxids=blast_results['subject_taxids'].str.split(';')).explode('subject_taxids')
    blast_results['subject_names'] = blast_results['subject_taxids'].apply(lambda x: type[x])
    blast_results['description'] = blast_results['subject_accession'].apply(lambda x: seq_description[x])
    
    blast_results.to_csv(args.outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-f','--taxonresult')
    parser.add_argument('-d','--descriptionfile')
    parser.add_argument('-b','--blastresult')
    parser.add_argument('-o','--outfile')
    args = parser.parse_args()
    main(args)