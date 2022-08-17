#!/usr/bin/env python3

import sys as sys
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd
import re
import argparse





#ref = "/scratch/sherrie.wang/hcv_round2/ref/HCV_ref_genomes_2019-8-5-19_coding_region_only_nospaces.fasta"
amplicons_1 = {'forward':'AGGTCTCGTAGACCGTGCATCATG','reverse':'CA[CTGA]GT[CTGA]AGGGTATCGATGAC'}
amplicons_2 = {'forward':'TATGA[CTGA]ACCCGCTG[CTGA]TTTGACTC','reverse':'GC[CTGA]GA[AGCT]TA[CTGA]CT[ACGT]GTCATAGCCTC'}
#samtools view A2-corens_alignment_filtered_sorted.bam | cut -f3 | head -n 1 > genotype_name.txt
#genotype_file = "genotype_name.txt"
#sample_name = "A2-corens"



def main(args):
    genotype_file = args.genotype_file
    sample_name = args.sample_name
    #read in genotype name
    with open(genotype_file) as f:
        genotype_name = f.readlines()

    genotype_name = genotype_name[0].strip()

    subname = genotype_name.split("|")[2] + '_' + genotype_name.split("|")[1]

    #read in ref HCV sequences

    with open(args.db, 'r') as input_file:
        ref_seqs = {}
        for line in input_file:
            if line[0] == '>' :
                header = line.strip().lstrip('>')
                header = header.split(' ',1)[0]
                ref_seqs[header] = ''
            else:
                ref_seqs[header] += line.strip()

    #forward sequence
    seq = ref_seqs[subname]

    #reverse sequence
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement_seq = "".join(complement.get(base, base) for base in seq)

    if args.ends == 'core':
        bed = []

        #region 1
        indx1 = re.search(amplicons_1['forward'],seq)

        reverse_amp = amplicons_1['reverse'][::-1]
        reverse_amp = reverse_amp.replace(']',';')
        reverse_amp = reverse_amp.replace('[',']')
        reverse_amp =reverse_amp.replace(';','[')

        indx2 = re.search(reverse_amp,complement_seq)


        if indx1 is None:
            start1 = 0
        else:
            start1 = indx1.span()[0]

        if indx2 is None:
            end1 = 403
        else:
            end1 = indx2.span()[1]


        bed.append([genotype_name, str(start1), str(end1),'5end','60','+'])

        with open(sample_name+'_amplicon.bed', 'w') as outf:
            outf.write('\t'.join(bed[0]) + '\n' )    

    #region 2
    if args.ends == 'ns5b':
        bed2 = []
        indx3 = re.search(amplicons_2['forward'],seq)

        reverse_amp = amplicons_2['reverse'][::-1]
        reverse_amp = reverse_amp.replace(']',';')
        reverse_amp = reverse_amp.replace('[',']')
        reverse_amp =reverse_amp.replace(';','[')

        indx4 = re.search(reverse_amp,complement_seq)

        if indx3 is None:
            #start2 = 7900
            start2 = 0
        else:
            start2 = indx3.span()[0]

        if indx4 is None:
            #end2 = 8500
            end2 = 388
        else:
            end2 = indx4.span()[1]


        bed2.append([genotype_name, str(start2), str(end2),'3end','60','+'])

        with open(sample_name+'_amplicon.bed', 'w') as outf:
        #outf.write('\t'.join(bed2[0]) +'\n' +'\t'.join(bed2[1]) + '\n' )
            outf.write('\t'.join(bed2[0])  + '\n' )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--db')
    parser.add_argument('--ends')
    parser.add_argument('-o', '--sample_name', help="output prefix")
    parser.add_argument('-g', '--genotype_file', help = "text file for genotype that is the same in the bam alignment file")
    args = parser.parse_args()
    main(args)
