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

    #corrected_reads1,corrected_reads2 = error_corrector(args['-o'], args['-f'], args['-r'])
    #normalize_reads1,normalize_reads2 = normalize_reads(args['-o'], args['-f'], args['-r'])
    contigs = assemble_contigs(args['-o'], args['-f'], args['-r'])
    #contigs = assemble_contigs(args['-o'], args['-f'], args['-r'])
    #contigs_cleaned = cut_amplicon_contigs(args['-o'],contigs)
    #contigs = assemble_contigs(args['-o'], args['-f'], args['-r'])
    blast_out = align_contigs_to_ref_seqs(args['-o'], contigs, args['-d'],1)
    blast_results = filter_alignments(args['-o'], blast_out, args['-c'], args['-i'])
    if args['-m'] == 'assemble':
        ref_seqs = write_best_contigs_fasta(args['-o'], blast_results, contigs)
        #best_ref = write_best_ref_seqs_fasta(args['-o'], blast_results, args['-d'])
    elif args['-m'] == 'align':
        ref_seqs = write_best_ref_seqs_fasta(args['-o'], blast_results, args['-d'])
    #bam_out = map_reads(args['-o'], ref_seqs, args['-f'], args['-r'])
    #vcf_out = call_variants(args['-o'], args['-D'], args['-q'], ref_seqs, bam_out)
    #consensus_seqs = make_consensus_seqs(args['-o'], bam_out, args['-D'], args['-c'], ref_seqs, vcf_out)
    #blast_out_consensus = align_contigs_to_ref_seqs(args['-o'],consensus_seqs,args['-d'],2)
    #do genotype calls
    #genotypes = genotype_call(args['-o'],blast_out_consensus,args['-c'], args['-i'])

    #print('\nWRITING REPORT...')
    #sequenced_bases = count_sequenced_bases_in_consensus_seqs(consensus_seqs)
    #consensus_seq_lengths = get_consensus_seq_lengths(consensus_seqs)
    #reads_mapped_to_consensus_seqs = count_reads_mapped_to_consensus_seqs(args['-o'], bam_out)
    #write_reports(args['-o'], sequenced_bases, consensus_seq_lengths, reads_mapped_to_consensus_seqs)
    #clean_headers(consensus_seqs)
    #if args['-g'] == 'yes':
    #    garbage_collection(args['-o'])
    #print('\nDONE: Analysis finished succesfully.\n')
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


def normalize_reads(output, fwd_reads,rev_reads):
    '''bbnorm.sh on the trimmed reads.
    Returns path to normalized reads file.'''
    print('normalize reads depth for assembly...') 
    normalize_reads1 = os.path.join(output, output + '_R1_normalized.fastq')
    normalize_reads2 = os.path.join(output, output + '_R2_normalized.fastq')

    #terminal_command = f'reformat.sh in1={fwd_reads} in2={rev_reads} out1={normalize_reads1} out2={normalize_reads2} samplerate=0.3'
    terminal_command = f'bbnorm.sh in={fwd_reads} in2={rev_reads} out={normalize_reads1} out2={normalize_reads2} target=100 min=5'
    error_msg = f'bbnorm.sh terminated with errors. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_bbnorm_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bbnorm_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    return normalize_reads1, normalize_reads2  

def error_corrector(output, fwd_reads, rev_reads):
    '''correct errors'''
    print('Error correcting reads')
    correct_out = os.path.join(output, output + '_spades_error_correct')
    terminal_command = f'spades.py --only-error-correction -1 {fwd_reads} -2 {rev_reads} -o {correct_out}'
    error_msg = f'spades terminated with errors while error correcting reads. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_spades_correct_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_spades_correct_stderr.txt')
    corrected_reads1 = os.path.join(correct_out, 'corrected',output + '_R1.trim.fastq.00.0_0.cor.fastq.gz' )
    corrected_reads2 = os.path.join(correct_out, 'corrected',output + '_R2.trim.fastq.00.0_0.cor.fastq.gz' )
    run(terminal_command, error_msg, stdout_file, stderr_file) 
    return corrected_reads1, corrected_reads2   


def assemble_contigs(output, fwd_reads, rev_reads, contig_type='scaffolds'):
    '''Assmebles contigs from fwd_reads and rev_reads FASTQ files. Sends output to spades_out_dir.
    Returns path to contigs FASTA file.'''
    print('Assembling reads into contigs...')
    spades_out = os.path.join(output, output + '_spades_results')
    terminal_command = f'spades.py --rnaviral --isolate -k 15 21 25 -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'
    #terminal_command = f'spades.py --careful -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'
    error_msg = f'spades terminated with errors while assembling reads into contigs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_spades_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_spades_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    old_contigs = os.path.join(spades_out, f'{contig_type}.fasta')
    if os.path.exists(old_contigs) == False:
        print(f'\nDONE: No {contig_type} assembled from reads.')
        old_contigs = os.path.join(spades_out, f'contigs.fasta')
    #    #exit(0)
    contigs = os.path.join(output, output + '_contigs.fa')
    sh.copy2(old_contigs, contigs)
    return contigs

def cut_amplicon_contigs(output, ref):
    ''''''
    print('cut adapters from contigs')
    cutadapt_out = os.path.join(output, output + '_contigs_cleaned.fa')
    #terminal_command = f'spades.py --rnaviral --isolate -1 {fwd_reads} -2 {rev_reads} -o {spades_out}'TATGAYACCCGCTGYTTTGACTC
    #terminal_command = f'cutadapt -g "AGGTCTCGTAGACCGTGCATCATG" --rc -o intermediate1.fa {contigs}'
    #terminal_command2 = f'cutadapt -b "GTCATCGATACCCTYACRTG" --rc -o intermediate2.fa intermediate1.fa'
    #terminal_command3 = f'cutadapt -b "CAYGTRAGGGTATCGATGAC" --rc -o intermediate3.fa intermediate2.fa'
    #terminal_command4 = f'cutadapt -b "GAGGCTATGACBAGRTAYTCNGC" --rc -o {cutadapt_out} intermediate3.fa'
    terminal_command = f'bbduk.sh in=${ref} out=test_primer_remove.fa literal=AGGTCTCGTAGACCGTGCATCATG,GTCATCGATACCCTYACRTG,CAYGTRAGGGTATCGATGAC,GAGGCTATGACBAGRTAYTCNGC rcomp=t ktrim=rl k=19 mink=9 hdist=1 tpe tbo'

    error_msg = f'cutadapt terminated with errors while assembling reads into contigs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_cutadapt_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_cutadapt_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    run(terminal_command2, error_msg, stdout_file, stderr_file)
    run(terminal_command3, error_msg, stdout_file, stderr_file)
    run(terminal_command4, error_msg, stdout_file, stderr_file)
    return cutadapt_out

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
    #blast_results['end'] = blast_results.apply(lambda row: row['sseqid'].split('|')[1], axis=1)
    # Discard alignments below minimum identity threshold

    blast_results['coverage'] = blast_results['qlen'] * 100 / blast_results['slen']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results['amplicon'] = blast_results['send'].apply(condition)
    #write out blast results to check all alignments above an identify threshold.
    blast_results.to_csv(os.path.join(output, output + '_blast_results_prefilter.csv'))
    # Keep only best alingments for each contig (by bitscore)
    blast_results = blast_results[(blast_results['qlen'] >= 300)]
    #blast_results = blast_results[(blast_results['qlen'] <= 425)] #increase this range to allow PCR runthrough products wiggle room
    #blast_results = blast_results[(blast_results['qlen'] <= 420) | (blast_results['qend'] - blast_results['qstart'] >= 300)]
    #blast_results = blast_results[blast_results['pident']>=min_id]
    best_bitscores = blast_results[['qseqid', 'bitscore']].groupby('qseqid').max().reset_index()
    blast_results = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore'])
    # Discard contigs whose best alignments are to multiple segments
    #segment_counts = blast_results[['qseqid', 'segment']].drop_duplicates()
    #segment_counts = segment_counts.groupby('qseqid').size().reset_index()
    #segment_counts = segment_counts[segment_counts[0]==1][['qseqid']]
    #blast_results = pd.merge(blast_results, segment_counts, on='qseqid')
    # Discard contigs whose best alignments are to multiple subtypes
    #subtype_counts = blast_results[['qseqid', 'subtype']].drop_duplicates()
    #subtype_counts = subtype_counts.groupby('qseqid').size().reset_index()
    #subtype_counts = subtype_counts[subtype_counts[0]==1][['qseqid']]
    #blast_results = pd.merge(blast_results, subtype_counts, on='qseqid')
    # Keep only alignments between contigs and ref seqs with median segment length
    #median_slen = blast_results[['qseqid', 'slen']].groupby('qseqid').quantile(0.5, interpolation='higher').reset_index()
    #blast_results = pd.merge(blast_results, median_slen, on=['qseqid', 'slen'])
    # Discard contigs that do not provide minimum coverage of a segment

    #blast_results = blast_results[blast_results['qlen'] * 100 / blast_results['slen'] >= min_cov]
    #blast_results = blast_results[blast_results['qlen'] * 100 / blast_results['slen'] <= 4.5]
    
    #best_bitscores = blast_results[['subtype', 'bitscore']].groupby('subtype').max().reset_index()
    #blast_results = pd.merge(blast_results, best_bitscores, on=['subtype', 'bitscore'])

    # De-duplicate sheet
    #cols = ['qseqid', 'sseqid', 'segment', 'subtype', 'slen','bitscore']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results = blast_results.drop_duplicates(['amplicon'], keep='first') #for each end, keep the contig with the best coverage
    
    #blast_results = blast_results[cols].drop_duplicates()
    #blast_results=blast_results.sort_values(by=['bitscore'],ascending =False)
    #seqids = list(blast_results['sseqid'].drop_duplicates())
    #if len(seqids) > 2: #if more than three subtype and segments are present, use the top 2
    #    blast_results = blast_results[blast_results['sseqid'].isin(seqids[0:2])]
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

    #blast_results.to_csv(os.path.join(output, output + '_filtered_blast_results.csv'))
    return blast_results


def write_best_contigs_fasta(output, blast_results, contigs):
    '''Looks up best contigs in contigs FASTA file and writes them to their own FASTA file.
    Returns path to best contigs FASTA file.'''
    # De-duplicate rows from contigs with best alignments to multiple ref seqs 
    cols = ['qseqid', 'segment','subtype', 'slen','qstart','qend','amplicon']
    
    blast_results = blast_results[cols].drop_duplicates(subset=['qseqid', 'segment','subtype', 'slen','amplicon'], keep='first')

    # Open contigs FASTA and load seqs into dict (key=seq header)
    with open(contigs, 'r') as input_file:
        contig_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>').rstrip(' rc')
                contig_seqs[header] = ''
            else:
                contig_seqs[header] += line.strip()

    # Create path for best contigs FASTA
    best_contigs = os.path.join(output, output + '_filtered_contigs.fa')
    # Write best contigs to FASTA
    with open(best_contigs, 'w') as output_file:
        contig_counter = 1
        for index in blast_results.index:
            contig_name = blast_results['qseqid'][index]
            segment = blast_results['segment'][index]
            subtype = blast_results['subtype'][index]
            segment_length = blast_results['slen'][index]
            amplicon = blast_results['amplicon'][index]
            header = f'>{output}_seq_{contig_counter}|{amplicon}|{segment}|{subtype}|{segment_length}'
            #if (len(contig_seqs[contig_name]) > 420) & (blast_results['qend'][index] - blast_results['qstart'][index] > 300):
            #    contig_seqs[contig_name] = contig_seqs[contig_name][(blast_results['qstart'][index]):(blast_results['qend'][index])]
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


def map_reads(output, ref_seqs, fwd_reads, rev_reads):
    '''Maps reads to ref seqs (either contigs or first alphabetical ref seq) with BWA mem. Filters, sorts, and indexes mappings with samtools.
    Returns path to filtered/sorted/indexed BAM file.'''
    print('Mapping reads...')
    terminal_command = f'bwa index {ref_seqs}'
    error_msg = f'bwa index terminated with errors. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_bwa_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bwa_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    sam_out = os.path.join(output, output + '_alignment.sam')
    terminal_command = f'bwa mem {ref_seqs} {fwd_reads} {rev_reads} > {sam_out}'
    error_msg = f'bwa mem terminated with errors while mapping reads. Please refer to /{output}/logs/ for output logs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_bwa_mem_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    bam_out = os.path.join(output, output + '_alignment_filtered_sorted.bam')
    terminal_command = f'samtools view -f 3 -F 2828 -q 30 -h {sam_out} | samtools sort -o {bam_out}'
    #terminal_command = f'samtools view -S -b {sam_out} | samtools sort -o {bam_out}'
    error_msg = f'samtools view/sort terminated with errors while filtering and sorting mapped reads. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)

    #for x in range(2):
    #    #check if bam is empty
    #    bam_empty_file = os.path.join(output, output + '_bam_empty.txt')
    #    terminal_command = f'samtools view {bam_out} | wc -l > {bam_empty_file}'
    #    error_msg = f'samtools view/sort terminated with errors while cat bam file. Please refer to /{output}/logs/ for output logs.'
    #    stdout_file = os.path.join(output, 'logs', output + '_samtools_view_stdout.txt')
    #    stderr_file = os.path.join(output, 'logs', output + '_samtools_view_stderr.txt')
    #    run(terminal_command, error_msg, stdout_file, stderr_file)
    #    with open(bam_empty_file) as f:
    #        bam_empty = f.readlines()
    #    bam_empty = int(bam_empty[0].strip())

    #    if bam_empty == 0: #if it is empty then remove the filter on proper pair
    #        bam_out = os.path.join(output, output + '_alignment_filtered_sorted.bam')
    #        terminal_command = f'samtools view -f 1 -F 2828 -q 30 -h {sam_out} | samtools sort -o {bam_out}'
    #        if x==1: #allow reads align to multiple segments due to high similarity, so ignore mapping quality 
    #            terminal_command = f'samtools view -f 1 -F 2828 -h {sam_out} | samtools sort -o {bam_out}'
    #        
    #        #terminal_command = f'samtools view -S -b {sam_out} | samtools sort -o {bam_out}'
    #        error_msg = f'samtools view/sort terminated with errors while filtering and sorting mapped reads. Please refer to /{output}/logs/ for output logs.'
    #        stdout_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stdout.txt')
    #        stderr_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stderr.txt')
    #        run(terminal_command, error_msg, stdout_file, stderr_file)

    terminal_command = f'samtools index {bam_out}'
    error_msg = f'samtools index terminated with errors while indexing mapped reads. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_samtools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_samtools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    return bam_out


def call_variants(output, min_depth, min_qual, ref_seqs, bam_out):
    '''Call variants with bcftools. Returns path to VCF file.'''
    print('Calling variants from mapped reads...')
    vcf_out = os.path.join(output, output + '_variants.vcf.gz')
    terminal_command = (f'bcftools mpileup -q {min_qual} -Q {min_qual}'
                        f' -m {min_depth} -Ou -f {ref_seqs} {bam_out}'
                        f' | bcftools call --ploidy 1 -M -mv -Oz -o {vcf_out}')
    error_msg = f'bcftools mpileup/call terminated with errors while calling variants. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_mpileup_call_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_mpileup_call_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    terminal_command = f'bcftools index {vcf_out}'
    error_msg = f'bcftools index terminated with errors while indexing variant calls. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    return vcf_out


def make_consensus_seqs(output, bam_out, min_depth, min_cov, ref_seqs, vcf_out):
    '''Apply variants to consensus seqs with bcftools. Returns path to consensus seqs FASTA file.'''
    print('Masking low coverage positions...')
    low_cov = os.path.join(output, output + '_low_cov.bed')
    terminal_command = f"bedtools genomecov -bga -ibam {bam_out} | awk '$4<{min_depth} {{print}}' > {low_cov}"
    error_msg = f'bedtools genomecov terminated with errors while masking low coverage positions. Please refer to /{output}/logs/ for output logs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_bedtools_genomecov_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    print('Generating consensus seqs...')
    consensus_seqs = os.path.join(output, output + '_consensus_seqs.fa')
    terminal_command = (f'cat {ref_seqs} | bcftools consensus -m {low_cov} {vcf_out}'
                        f' | seqtk seq -l 0 > {consensus_seqs}')
    error_msg = f'bcftools consensus terminated with errors while applying variants. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_bcftools_index_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_bcftools_index_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    #terminal_command = f"cat {ref_seqs} | grep '>' | cut -d'_' -f1 | uniq | sed 's/>//g' > segment_order"
    #sp.run(terminal_command,shell=True)
    #sp.run(f"(echo order ; cat segment_order) > segment_order.csv", shell=True)
    #segment_order = pd.read_csv("segment_order.csv")
    #segment_order = segment_order['order'].tolist()
    #segment_order = ['1a', '1b', '1c', '1d', '1e', '1g', '1h', '1j','1k','1m','1n','1o','1','2a','2b','2c','2e','2f','2i','2j','2k','2m','2q','2r','2t','2u','2','3a','3b','3d','3e','3g','3h','3i','3k','3','4a','4b','4c','4d','4f','4g','4k','4l','4m','4n','4o','4p','4q','4r','4s','4t','4v','']
    #header_order = sorted(seqs.keys(), key=lambda s: segment_order.index(s.split('|')[2]))
    #seqs = {header: seqs[header] for header in header_order}
    # Remove seqs that where the number of sequenced bases does not exceed the min coverage of the median ref seq length
    ref_seq_length = lambda header: int(header.split('|')[-1])
    sequenced_bases = lambda seq: len([base for base in seq if base in 'ATGC'])
    sufficient_cov = lambda header, seq: True if sequenced_bases(seq) * 100 / ref_seq_length(header) >= min_cov else False
    seqs = {header: seq for header, seq in seqs.items() if sufficient_cov(header, seq) == True}
    # Check if any consensus seqs remain
    if len(seqs) == 0:
        print('\nDONE: No consensus seqs exceeding minimum coverage of segment.\n')
        exit(0)
    # Write out consensus seqs
    with open(consensus_seqs, 'w') as output_file:
        for header, seq in seqs.items():
            output_file.write(header + '\n')
            output_file.write(seq + '\n')
    return consensus_seqs


def count_sequenced_bases_in_consensus_seqs(consensus_seqs):
    print('Counting sequenced bases in consensus seqs...')
    check_input_exists([consensus_seqs])
    check_input_empty([consensus_seqs])
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    sequenced_bases = {header: len(seq) - seq.count('N') for header, seq in seqs.items()}
    sequenced_bases = pd.DataFrame.from_dict(sequenced_bases, orient='index').reset_index()
    sequenced_bases.columns = ['consensus_seq', 'sequenced_bases']
    return sequenced_bases


def get_consensus_seq_lengths(consensus_seqs):
    print('Getting lengths of consensus seqs...')
    check_input_exists([consensus_seqs])
    check_input_empty([consensus_seqs])
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    consensus_seq_lengths = {header: len(seq) for header, seq in seqs.items()}
    consensus_seq_lengths = pd.DataFrame.from_dict(consensus_seq_lengths, orient='index').reset_index()
    consensus_seq_lengths.columns = ['consensus_seq', 'seq_length']
    return consensus_seq_lengths


def count_reads_mapped_to_consensus_seqs(output, bam_out):
    idxstats = os.path.join(output, output + '_reads_mapped_to_consensus_seqs.tsv')
    terminal_command = f'samtools idxstats {bam_out} > {idxstats}'
    error_msg = f'samtools idxstats terminated with errors while counting reads mapped to consensus seqs. Please refer to /{output}/logs/ for output logs.'
    stdout_file = None
    stderr_file = os.path.join(output, 'logs', output + '_samtools_idxstats_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)
    cols = ['consensus_seq', 'seq_length', 'mapped_reads', 'unmapped_reads']
    reads_mapped_to_consensus_seqs = pd.read_csv(idxstats, sep='\t', names=cols).replace('*', np.nan).dropna()
    cols = ['consensus_seq', 'mapped_reads']
    reads_mapped_to_consensus_seqs = reads_mapped_to_consensus_seqs[cols]
    return reads_mapped_to_consensus_seqs

def genotype_call(output, consensus_blast_out, min_cov, min_id):
    '''Find best contig for each genome segment. Returns datasheet with best contigs.'''
    print('Filtering alignments...')
    cols = 'qseqid sseqid pident qlen slen mismatch gapopen qstart qend sstart send bitscore'.split(' ')
    blast_results = pd.read_csv(consensus_blast_out, sep='\t', names=cols)
    blast_results['segment'] = blast_results.apply(lambda row: row['sseqid'].split('_')[1], axis=1)
    blast_results['subtype'] = blast_results.apply(lambda row: row['sseqid'].split('_')[0], axis=1)
    blast_results = blast_results[blast_results['pident'] >= 90]
    #blast_results['end'] = blast_results.apply(lambda row: row['sseqid'].split('|')[1], axis=1)
    # Discard alignments below minimum identity threshold

    blast_results['coverage'] = blast_results['qlen'] * 100 / blast_results['slen']
    blast_results = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results['amplicon'] = blast_results.apply(lambda row: row['qseqid'].split('|')[1], axis=1)
    #blast_results = blast_results[blast_results['pident']>=min_id]

    best_bitscores = blast_results[['qseqid','bitscore']].groupby(['qseqid']).head(10).reset_index()
    blast_results1 = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore'])
    #blast_results1 = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results1.to_csv(os.path.join(output, output + '_genotype_calls.csv'))

    return blast_results1

def write_reports(output, sequenced_bases, consensus_seq_lengths, reads_mapped_to_consensus_seqs):
    print('Compiling data for report...')
    print(reads_mapped_to_consensus_seqs)
    print(sequenced_bases)
    print(consensus_seq_lengths)
    report = pd.merge(reads_mapped_to_consensus_seqs, sequenced_bases, on='consensus_seq')
    report = pd.merge(report, consensus_seq_lengths, on='consensus_seq')
    report['amplicon'] = report.apply(lambda row: row['consensus_seq'].split('|')[1], axis=1)
    report['segment'] = report.apply(lambda row: row['consensus_seq'].split('|')[2], axis=1)
    report['subtype'] = report.apply(lambda row: row['consensus_seq'].split('|')[3], axis=1)
    report['ref_seq_length'] = report.apply(lambda row: int(row['consensus_seq'].split('|')[-1]), axis=1)
    report['consensus_seq'] = report.apply(lambda row: row['consensus_seq'].split('|')[0], axis=1)
    report['segment_cov'] = round(report['sequenced_bases'] * 100 / report['ref_seq_length'], 2)
    print(report)
    print('Writing consensus seq report...')
    cols = ['consensus_seq','amplicon', 'segment', 'subtype', 'mapped_reads', 'seq_length', 'sequenced_bases', 'segment_cov']
    contig_report = report[cols].drop_duplicates()
    report_path = os.path.join(output, output + '_consensus_seqs_report.tsv')
    contig_report.to_csv(report_path, sep='\t', index=False)


def clean_headers(consensus_seqs):
    with open(consensus_seqs, 'r') as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip()
                seqs[header] = ''
            else:
                seqs[header] += line.strip()
    with open(consensus_seqs, 'w') as output_file:
        for header, seq in seqs.items():
            header = '|'.join(header.split('|')[:-1])
            output_file.write(header + '\n')
            output_file.write(seq + '\n')


def garbage_collection(output):
    spades_out = os.path.join(output, output + '_spades_results')
    sh.rmtree(spades_out)
    blast_out = os.path.join(output, output + '_contigs_blast_results.tsv')
    os.remove(blast_out)
    #contigs = os.path.join(output, output + '_contigs.fa')
    #os.remove(contigs)
    sam_out = os.path.join(output, output + '_alignment.sam')
    os.remove(sam_out)
    ref_seqs = os.path.join(output, output + '_ref_seqs_for_mapping.fa')
    vcf_out = os.path.join(output, output + '_variants.vcf.gz')
    low_cov = os.path.join(output, output + '_low_cov.bed')
    files = [ref_seqs + suffix for suffix in ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']]
    files += [vcf_out + suffix for suffix in ['', '.csi']]
    files += [low_cov]
    for file in files:
        os.remove(file)
    bam_out = os.path.join(output, output + '_alignment_filtered_sorted.bam')
    idxstats = os.path.join(output, output + '_reads_mapped_to_consensus_seqs.tsv')
    files = [bam_out + suffix for suffix in ['.bai']] + [idxstats]
    for file in files:
        os.remove(file)


if __name__ == '__main__':
    main()
