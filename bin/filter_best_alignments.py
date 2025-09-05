#!/usr/bin/env python3

import sys as sys
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd


def main():
    args = parse_args(sys.argv)

    blast_results = filter_alignments(args['-o'], args['-f'], args['-c'], args['-i'])
    if args['-m'] == 'assemble':
        ref_seqs = write_best_contigs_fasta(args['-o'], blast_results, args['-t'])
    elif args['-m'] == 'align':
        ref_seqs = write_best_ref_seqs_fasta(args['-o'], blast_results, args['-d'])

    exit(0)


def parse_args(args):
    arg_values = {}
    for arg_1, arg_2 in zip(args[1:-1], args[2:]):
        if arg_1[0] == '-':
            if arg_2[0] != '-':
                arg_values[arg_1] = arg_2
            else:
                arg_values[arg_1] = ''
    if args[-1][0] == '-':
        arg_values[args[-1]] = ''
    # Set defaults, mins, and maxs
    required_args = {'-f','-d', '-m', '-o','-t'}
    arg_value_types = {'-f': str, '-d': str, '-m': str, '-o': str, '-c': float, '-i': float, '-t':str}
    min_arg_values = { '-c': 0, '-i': 0}
    max_arg_values = {'-c': 100, '-i': 100}
    default_arg_values = { '-c': 1, '-i': 90}
    # Check if all required arguments were provided
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys() or arg_values[required_arg] == '':
            missing_args = missing_args | {required_arg}
    if missing_args != set():
        print(f'\nERROR: Values must be provided for the argument following arguments: {", ".join(sorted(missing_args))}')
        exit(1)
    # Check if unrecognized arguments were provided
    recognized_args = required_args | set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args = unrecognized_args | {provided_arg}
    if unrecognized_args != set():
        print(f'\nERROR: The following arguments are not recognized: {", ".join(sorted(unrecognized_args))}')
        exit(1)
    # Check if arguments were provided without values
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args = empty_args | {arg}
    if empty_args != set():
        print(f'\nERROR: The following arguments were provided without values: {", ".join(sorted(empty_args))}')

        exit(1)
    # Check if provided values are of the correct type
    for arg, value in arg_values.items():
        try:
            arg_values[arg] = arg_value_types[arg](value)
        except ValueError:
            print(f'\nERROR: Value for argument {arg} must be of type {str(arg_value_types[arg].__name__)}')

            exit(1)
    # Check if provided values are within the correct range
    for arg, value in arg_values.items():
        if arg in min_arg_values.keys() and value <= min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be greater than {min_arg_values[arg]}')
    
            exit(1)
        if arg in max_arg_values.keys() and value >= max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be less than {max_arg_values[arg]}')

            exit(1)        
    # Assign default values to unspecified arguments
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value 
    # Return keyword args and their values
    return arg_values
         

def condition(x):
    if x>7500:
        return "ns5b"
    elif x<=850:
        return "core"
    else:
        return 'weird'

def filter_alignments(output, blast_out, min_cov, min_id):
    '''Find best contig for each genome segment. Returns datasheet with best contigs.'''
    print('Filtering alignments...')
    cols = 'qseqid sseqid pident qlen slen mismatch gapopen qstart qend sstart send bitscore'.split(' ')
    blast_results = pd.read_csv(blast_out, sep='\t', names=cols)
    # Annotate alignments with segment and subtype
    blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('_')[1], axis=1)
    blast_results['subtype'] = blast_results.apply(lambda row: row['sseqid'].split('_')[0], axis=1)
    # Discard alignments below minimum identity threshold

    blast_results['coverage'] = blast_results['qlen'] * 100 / blast_results['slen']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results['amplicon'] = blast_results['send'].apply(condition)
    #write out blast results to check all alignments above an identify threshold.
    blast_results.to_csv(os.path.join(output + '_blast_results_prefilter.csv'))
    # Keep only best alingments for each contig (by bitscore)
    blast_results = blast_results[(blast_results['qlen'] >= 300)]
    best_bitscores = blast_results[['qseqid', 'bitscore','amplicon']].groupby(['qseqid','amplicon']).max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore','amplicon'])
    best_bitscores = blast_results[['amplicon','subtype', 'bitscore']].groupby(['amplicon','subtype']).max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['amplicon','subtype', 'bitscore'])
    #best_bitscores = blast_results.groupby(['amplicon','subtype']).head(1).reset_index()

    #only allow one contig to one subtype, with the same bitscore
    cols=['amplicon','subtype','bitscore']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results = blast_results.drop_duplicates(cols, keep='first')

    # De-duplicate sheet
    cols = ['qseqid', 'subtype', 'amplicon']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results = blast_results.drop_duplicates(cols, keep='first') #for each end, keep the contig with the best coverage
    blast_results = blast_results[(abs(blast_results['qend'] - blast_results['qstart']) >= 200)]

    print(blast_results)
    blast_results.to_csv(os.path.join(output + '_filtered_blast_results.csv'))

    #parse blast results
    counts = blast_results['subtype'].value_counts().reset_index()
    print(counts)
    counts=counts.rename(columns = {'subtype' : 'counts','index':'subtype'})
    total=counts[['counts']].sum()

    counts['prop'] = counts['counts'].apply(lambda x: x/total)
    max_bitscore = blast_results[['subtype','bitscore']].groupby('subtype').max().reset_index()
    bit_df = pd.merge(counts,max_bitscore,how ='left', on='subtype')
    bit_df=bit_df.rename(columns = {'bitscore': 'max_bitscore'})
    bit_df.to_csv(os.path.join(output + '_max_bitscores_per_subtype.csv'))

    #blast_results = blast_results.head(5)
    print(blast_results)
    if len(blast_results) == 0:
        print('DONE: No valid contigs found.')
        exit(0)

    return blast_results


def write_best_contigs_fasta(output, blast_results, contigs):
    '''Looks up best contigs in contigs FASTA file and writes them to their own FASTA file.
    Returns path to best contigs FASTA file.'''
    # De-duplicate rows from contigs with best alignments to multiple ref seqs 
    #TODO: add a condition to not write the contig if the same sequence already 
    cols = ['qseqid', 'segment','subtype', 'slen','qstart','qend','amplicon']
    
    blast_results = blast_results[cols].drop_duplicates(subset=['qseqid', 'segment','subtype', 'slen','amplicon'], keep='first')

    # Open contigs FASTA and load seqs into dict (key=seq header)
    with open(contigs, 'r') as input_file:
        contig_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>').rstrip(' rc').split(' ')[0]
                contig_seqs[header] = ''
            else:
                contig_seqs[header] += line.strip()

    # Create path for best contigs FASTA
    best_contigs = os.path.join(output + '_filtered_contigs.fa')
    # Write best contigs to FASTA
    with open(best_contigs, 'w') as output_file:
        contig_counter = 1
        seqs_to_write = {}
        for index in blast_results.index:
            contig_name = blast_results['qseqid'][index]
            segment = blast_results['segment'][index]
            subtype = blast_results['subtype'][index]
            segment_length = blast_results['slen'][index]
            amplicon = blast_results['amplicon'][index]
            header = f'>{output}_seq_{contig_counter}|{amplicon}|{segment}|{subtype}|{segment_length}'
            if contig_seqs[contig_name] in seqs_to_write.values() and len(contig_seqs[contig_name]) < 600 :
                continue
            
        
            seqs_to_write[contig_name] = ''
            seqs_to_write[contig_name] += contig_seqs[contig_name]
            output_file.write(header + '\n')
            output_file.write(contig_seqs[contig_name] + '\n')

            contig_counter += 1
    return best_contigs


def write_best_ref_seqs_fasta(output, blast_results, ref_seqs_db):
    '''Looks up best ref seqs in ref seqs DB FASTA file and writes them to their own FASTA file.
    Returns path to best ref seqs FASTA file.'''
    # Check if multiple HA or NA subtypes are present
    #HA_subtypes = blast_results[blast_results['segment']=='HA']['subtype'].unique()
    #NA_subtypes = blast_results[blast_results['segment']=='NA']['subtype'].unique()
    #if len(HA_subtypes) > 1 or len(NA_subtypes) > 1:
        #print('WARNING: Multiple HA or NA subtypes detected. Internal segment sequences are not generated for mixed infections in align mode.')
        #blast_results = blast_results[blast_results['segment'].isin(['HA', 'NA'])]
    # Choose ref seqs with max bitscore for each segment/subtype combination
    best_bitscores = blast_results[['segment', 'subtype', 'bitscore','amplicon']].groupby(['segment', 'subtype','amplicon']).max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['segment', 'subtype', 'bitscore','amplicon'])
    # Chose ref seqs with median length for each segment/subtype combination
    median_lengths = blast_results[['segment', 'subtype', 'slen','amplicon']].groupby(['segment', 'subtype','amplicon']).quantile(0.5, interpolation='higher').reset_index()
    blast_results = pd.merge(blast_results, median_lengths, on=['segment', 'subtype', 'slen','amplicon'])
    # Choose first alphabetical ref seq for each segment/subtype combination
    first_ref_seqs = blast_results[['sseqid', 'segment', 'subtype','amplicon']].groupby(['segment', 'subtype','amplicon']).min().reset_index()
    blast_results = pd.merge(blast_results, first_ref_seqs, on=['sseqid', 'segment', 'subtype','amplicon'])
    # De-duplicate alignments
    cols = ['sseqid', 'segment', 'subtype', 'slen','amplicon']
    blast_results = blast_results[cols].drop_duplicates()
    # Open ref seqs DB FASTA and load seqs into dict (key=seq header)
    with open(ref_seqs_db, 'r') as input_file:
        ref_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                header = header.split(' ',1)[0]
                ref_seqs[header] = ''
            else:
                ref_seqs[header] += line.strip()
    # Create path for best ref seqs FASTA  
    best_ref_seqs = os.path.join(output + '_best_ref.fa')
    # Write best ref seqs to FASTA
    with open(best_ref_seqs, 'w') as output_file:
        ref_seq_counter = 1
        for index in blast_results.index:
            ref_seq_name = blast_results['sseqid'][index]
            segment = blast_results['segment'][index]
            subtype = blast_results['subtype'][index]
            segment_length = blast_results['slen'][index]
            header = f'>{output}_seq_{ref_seq_counter}|{segment}|{subtype}|{segment_length}'
            output_file.write(header + '\n')
            output_file.write(ref_seqs[ref_seq_name] + '\n')
            ref_seq_counter += 1
    return best_ref_seqs

if __name__ == '__main__':
    main()
