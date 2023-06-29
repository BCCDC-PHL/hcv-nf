#!/usr/bin/env python3

import sys as sys
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd


def main():
    version = '0.0.2'
    args = parse_args(sys.argv, version)
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer\n')
    if args['-m'] not in ['assemble', 'align']:
        print(f'\nERROR: "{args["-m"]}" is not a valid FluViewer mode. Must choose "assemble" or "align".\n')
        exit(1)
    print(f'Analysis mode: {args["-m"]}')
    print(f'Ref seq DB: {args["-d"]}')
    print()
    print(f'Analyzing library {args["-o"]}...')
    check_input_exists([args['-f'], args['-r'], args['-d']])
    check_input_empty([args['-f'], args['-r'], args['-d']])
    print(f'Fwd reads: {args["-f"]}')
    print(f'Rev reads: {args["-r"]}')
    print()
    make_out_dir(args['-o'])
    print('\nGENERATING CONSENSUS SEQS...')

    contigs = assemble_contigs(args['-o'],args['-f'], args['-r'])

    blast_out = align_contigs_to_ref_seqs(args['-o'], contigs, args['-d'],1)

    blast_results = filter_alignments(args['-o'], blast_out, args['-c'], args['-i'])
    if args['-m'] == 'assemble':
        ref_seqs = write_best_contigs_fasta(args['-o'], blast_results, contigs)
    elif args['-m'] == 'align':
        ref_seqs = write_best_ref_seqs_fasta(args['-o'], blast_results, args['-d'])

    exit(0)


def parse_args(args, version):
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
    required_args = {'-f', '-r', '-d', '-m', '-o'}
    arg_value_types = {'-f': str, '-r': str, '-d': str, '-m': str, '-o': str, '-D': int, '-q': int, '-c': float, '-i': float, '-g': str,'--adapter_sequence':str, '--adapter_sequence_2':str}
    min_arg_values = {'-D': 1, '-q': 0, '-c': 0, '-i': 0}
    max_arg_values = {'-c': 100, '-i': 100}
    default_arg_values = {'-D': 20, '-q': 30, '-c': 1, '-i': 90, '-g': 'yes','--adapter_sequence':'', '--adapter_sequence_2':''}
    # Check if all required arguments were provided
    missing_args = set()
    for required_arg in required_args:
        if required_arg not in arg_values.keys() or arg_values[required_arg] == '':
            missing_args = missing_args | {required_arg}
    if missing_args != set():
        print(f'\nERROR: Values must be provided for the argument following arguments: {", ".join(sorted(missing_args))}')
        print_usage(version)
        exit(1)
    # Check if unrecognized arguments were provided
    recognized_args = required_args | set(arg_value_types.keys())
    unrecognized_args = set()
    for provided_arg in arg_values.keys():
        if provided_arg not in recognized_args:
            unrecognized_args = unrecognized_args | {provided_arg}
    if unrecognized_args != set():
        print(f'\nERROR: The following arguments are not recognized: {", ".join(sorted(unrecognized_args))}')
        print_usage(version)
        exit(1)
    # Check if arguments were provided without values
    empty_args = set()
    for arg, value in arg_values.items():
        if value == '':
            empty_args = empty_args | {arg}
    if empty_args != set():
        print(f'\nERROR: The following arguments were provided without values: {", ".join(sorted(empty_args))}')
        print_usage(version)
        exit(1)
    # Check if provided values are of the correct type
    for arg, value in arg_values.items():
        try:
            arg_values[arg] = arg_value_types[arg](value)
        except ValueError:
            print(f'\nERROR: Value for argument {arg} must be of type {str(arg_value_types[arg].__name__)}')
            print_usage(version)
            exit(1)
    # Check if provided values are within the correct range
    for arg, value in arg_values.items():
        if arg in min_arg_values.keys() and value <= min_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be greater than {min_arg_values[arg]}')
            print_usage(version)
            exit(1)
        if arg in max_arg_values.keys() and value >= max_arg_values[arg]:
            print(f'\nERROR: Value for argument {arg} must be less than {max_arg_values[arg]}')
            print_usage(version)
            exit(1)        
    # Assign default values to unspecified arguments
    for arg, value in default_arg_values.items():
        if arg not in arg_values.keys():
            arg_values[arg] = value
    # Return keyword args and their values
    return arg_values


def print_usage(version):
    print(f'\nFluViewer v{version}')
    print('https://github.com/KevinKuchinski/FluViewer\n')
    print('Usage: FluViewer -f <path_to_fwd_reads> -r <path_to_rev_reads> -d <path_to_db_file> -o <output_name> -m <mode> [-D <min_depth> -q <min_qual> -c <min_cov> -i <min_id> -g <garbage_collection>]\n')
    print('Required arguments:')
    print(' -f : path to FASTQ file containing forward reads')
    print(' -r : path to FASTQ file containing reverse reads')
    print(' -d : path to FASTA file containing FluViewer database')
    print(' -o : output name')
    print(' -m : run mode ("align" or "assemble")\n')
    print('Optional arguments:')
    print(' -D : minimum read depth for base calling (default = 20)')
    print(' -q : Minimum PHRED score for base quality and mapping quality (default = 30)')
    print(' -c : Minimum coverage of database reference sequence by contig (percentage, default = 25)')
    print(' -i : Minimum nucleotide sequence identity between database reference sequence and contig (percentage, default = 90)')
    print(' -g : Collect garbage and remove intermediate files ("yes" or "no", default = "yes")\n')


def run(terminal_command, error_msg, stdout_file, stderr_file):
    '''Runs terminal command and directs stdout and stderr into indicated files.'''
    log_files = {}
    for file, dest in zip([stdout_file, stderr_file], ['stdout', 'stderr']):
        if file != None:
            log_files[dest] = open(file, 'w')
            log_files[dest].write('*' * 80 + '\n')
            log_files[dest].write('Terminal command:\n')
            log_files[dest].write(terminal_command + '\n')
            log_files[dest].write('*' * 80 + '\n')
        else:
            log_files[dest] = None
    completed_process = sp.run(terminal_command, stdout=log_files['stdout'], stderr=log_files['stderr'], shell=True)        
    for file in zip([stdout_file, stderr_file], ['stdout', 'stderr']):
        if file != None:
            log_files[dest].close()
    if completed_process.returncode != 0:
        print('\nERROR:', error_msg)
        exit(1)
    return completed_process


def check_input_exists(file_list):
    '''Checks if input files exist.'''
    for file in file_list:
        if os.path.exists(file) == False:
            print('\nERROR: Input file does not exist:')
            print(file)
            exit(1)

def check_input_empty(file_list):
    '''Checks if input files are not empty.'''
    for file in file_list:
        if os.path.getsize(file) == 0:
            print('\nERROR: Input file is empty:')
            print(file)
            exit(1)            


def make_out_dir(out_dir):
    '''Creates output dir and a dir for logs within.'''
    print('Creating directory for output...')
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)
        logs_path = os.path.join(out_dir, 'logs')
        os.mkdir(logs_path)
    else:
        if os.path.isdir(out_dir) == False:
            print('\nERROR: Cannot create output directory because a file with that name already exists.')
            exit(1)
        if not os.listdir(out_dir):
            print('\nWARNING: Output directory already exists but is empty. Analysis will continue.')
            os.mkdir(out_dir)
            logs_path = os.path.join(out_dir, 'logs')
            os.mkdir(logs_path)
        else:
            print('\nERROR: Output directory already exists and is not empty.')
            exit(1)


def assemble_contigs(output, fwd_reads, rev_reads, contig_type='scaffolds'):
    '''Assmebles contigs from fwd_reads and rev_reads FASTQ files. Sends output to spades_out_dir.
    Returns path to contigs FASTA file.'''
    print('Assembling reads into contigs...')

    spades_out = os.path.join(output, output + '_spades_results')
    terminal_command = f'spades.py --rnaviral --isolate -k 15 21 25 -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'
    #terminal_command = f'spades.py -k 15,21,25 --careful --only-assembler -1 {sampled1} -2 {sampled2} -o {spades_out}'
    error_msg = f'spades terminated with errors while assembling reads into contigs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_spades_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_spades_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    old_contigs_1 = os.path.join(spades_out, f'{contig_type}.fasta')
    if os.path.exists(old_contigs_1) == False:
        print(f'\nDONE: No {contig_type} assembled from reads.')
        old_contigs_1 = os.path.join(spades_out, f'contigs.fasta')
    #    #exit(0)
    contigs = os.path.join(output, output + '_contigs.fa')
    #sh.copy2(old_contigs, contigs)
    sampled1 = os.path.join(output + '_sampled_R1.fq ')
    sampled2 = os.path.join(output + '_sampled_R2.fq ')
    terminal_command = (f'reformat.sh in1={fwd_reads} in2={rev_reads} out1={sampled1} out2={sampled2} samplerate=0.1')
    #error_msg = f'subsampling with reformat terminated with errors Please refer to /{output}/logs/ for output logs.'
    #stdout_file = os.path.join(output, 'logs', output + '_subsampling_stdout.txt')
    #stderr_file = os.path.join(output, 'logs', output + '_subsampling_stderr.txt')

    completed_process = sp.run(terminal_command, shell=True)
    if completed_process.returncode != 0: #if subsampling failed then run spades careful without subsampling 
        spades_out = os.path.join(output, output + '_spades_results_careful')
        terminal_command = f'spades.py -k 15,21,25 --careful --only-assembler -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'
        completed_process = sp.run(terminal_command, shell=True)
        if completed_process.returncode != 0: #if spades careful failed, then return old_contigs_1
            sh.copy2(old_contigs_1, contigs)
            return contigs

    #run(terminal_command, error_msg, stdout_file, stderr_file)
    spades_out = os.path.join(output, output + '_spades_results_careful')
    #otherwise if subsampling is successful, then run spades with careful with subsampled data
    terminal_command = f'spades.py -k 15,21,25 --careful --only-assembler -1 {sampled1} -2 {sampled2} -o {spades_out}'
    error_msg = f'spades terminated with errors while assembling reads into contigs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_spades_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_spades_stderr.txt')

    log_files = {}
    for file, dest in zip([stdout_file, stderr_file], ['stdout', 'stderr']):
        if file != None:
            log_files[dest] = open(file, 'w')
            log_files[dest].write('*' * 80 + '\n')
            log_files[dest].write('Terminal command:\n')
            log_files[dest].write(terminal_command + '\n')
            log_files[dest].write('*' * 80 + '\n')
        else:
            log_files[dest] = None
    completed_process = sp.run(terminal_command, stdout=log_files['stdout'], stderr=log_files['stderr'], shell=True)  
    #completed_process = sp.run(terminal_command, shell=True)
        #run(terminal_command, error_msg, stdout_file, stderr_file)
    if completed_process.returncode != 0: # if spades with careful failed then return old_contigs_1
        sh.copy2(old_contigs_1, contigs)
        return contigs
    
    old_contigs_2 = os.path.join(spades_out, f'{contig_type}.fasta')
    if os.path.exists(old_contigs_2) == False:
        print(f'\nDONE: No {contig_type} assembled from reads.')
        old_contigs_2 = os.path.join(spades_out, f'contigs.fasta')
    #    #exit(0)
    contigs = os.path.join(output, output + '_contigs.fa')
    #sh.copy2(old_contigs, contigs)
    terminal_command = f'cat {old_contigs_1} {old_contigs_2} > {contigs}'
    error_msg = f'cat terminated with errors while concat contigs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_cat_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_cat_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)

    return contigs


def align_contigs_to_ref_seqs(output, contigs, ref_seqs_db,time):
    '''Align contigs to reference sequences with BLASTn. Returns path to BLASTn results in TSV file.'''
    print(f'Aligning contigs to ref seqs in {ref_seqs_db}...')
    if any([os.path.exists(ref_seqs_db + '.' + suffix) == False for suffix in ['nhr', 'nin' , 'nsq']]):
        print(f'WARNING: blastn db files do not exist for {ref_seqs_db}. Creating blastn db files...')
        terminal_command = (f'makeblastdb -in {ref_seqs_db} -dbtype nucl')
        error_msg = f'blastn terminated with errors while making db for ref seqs. Please refer to /{output}/logs/ for output logs.'
        stdout_file = os.path.join(output, 'logs', output + '_make_blast_db_stdout.txt')
        stderr_file = os.path.join(output, 'logs', output + '_make_blast_db_stderr.txt')
        run(terminal_command, error_msg, stdout_file, stderr_file)
    if time == 1:
        blast_out = os.path.join(output, output + '_contigs_blast_results.tsv')
    if time == 2:
        blast_out = os.path.join(output, output + '_consensus_blast_results.tsv')
    terminal_command = (f'blastn -query {contigs} -db {ref_seqs_db} -outfmt '
                        f' "6 qseqid sseqid pident qlen slen mismatch gapopen qstart qend sstart send bitscore" > {blast_out}')
    error_msg = f'blastn terminated with errors while aligning contigs to ref seqs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_contigs_blast_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    if os.path.getsize(blast_out) == 0:
        print('\nDONE: No contigs aligned to ref seqs.')
        exit(0)        
    return blast_out

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
    blast_results.to_csv(os.path.join(output, output + '_blast_results_prefilter.csv'))
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
    blast_results.to_csv(os.path.join(output, output + '_filtered_blast_results.csv'))

    #parse blast results
    counts = blast_results['subtype'].value_counts().reset_index()
    counts=counts.rename(columns = {'subtype' : 'counts','index':'subtype'})
    total=counts[['counts']].sum()
    counts['prop'] = counts['counts'].apply(lambda x: x/total)
    max_bitscore = blast_results[['subtype','bitscore']].groupby('subtype').max().reset_index()
    bit_df = pd.merge(counts,max_bitscore,how ='left', on='subtype')
    bit_df=bit_df.rename(columns = {'bitscore': 'max_bitscore'})
    bit_df.to_csv(os.path.join(output, output + '_max_bitscores_per_subtype.csv'))

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
    best_contigs = os.path.join(output, output + '_filtered_contigs.fa')
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
    best_bitscores = blast_results[['segment', 'subtype', 'bitscore']].groupby(['segment', 'subtype']).max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['segment', 'subtype', 'bitscore'])
    # Chose ref seqs with median length for each segment/subtype combination
    median_lengths = blast_results[['segment', 'subtype', 'slen']].groupby(['segment', 'subtype']).quantile(0.5, interpolation='higher').reset_index()
    blast_results = pd.merge(blast_results, median_lengths, on=['segment', 'subtype', 'slen'])
    # Choose first alphabetical ref seq for each segment/subtype combination
    first_ref_seqs = blast_results[['sseqid', 'segment', 'subtype']].groupby(['segment', 'subtype']).min().reset_index()
    blast_results = pd.merge(blast_results, first_ref_seqs, on=['sseqid', 'segment', 'subtype'])
    # De-duplicate alignments
    cols = ['sseqid', 'segment', 'subtype', 'slen']
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
    best_ref_seqs = os.path.join(output, output + '_best_ref.fa')
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
