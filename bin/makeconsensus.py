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

    bam_out = map_reads(args['-o'], args['-s'], args['-f'], args['-r'])
    vcf_out = call_variants(args['-o'], args['-D'], args['-q'], args['-s'], bam_out)
    consensus_seqs = make_consensus_seqs(args['-o'], bam_out, args['-D'], args['-c'], args['-s'], vcf_out)
    
    #do genotype calls
    consensus_seqs_no_primers = cut_primers(args['-o'],consensus_seqs)

    #blast_out_consensus = align_contigs_to_ref_seqs(args['-o'],consensus_seqs_no_primers,args['-d'],2)
    #genotypes = genotype_call(args['-o'],blast_out_consensus,args['-c'], args['-i'])

    print('\nWRITING REPORT...')
    sequenced_bases = count_sequenced_bases_in_consensus_seqs(consensus_seqs)
    consensus_seq_lengths = get_consensus_seq_lengths(consensus_seqs)
    reads_mapped_to_consensus_seqs = count_reads_mapped_to_consensus_seqs(args['-o'], bam_out)
    write_reports(args['-o'], sequenced_bases, consensus_seq_lengths, reads_mapped_to_consensus_seqs)
    clean_headers(consensus_seqs)
    if args['-g'] == 'yes':
        garbage_collection(args['-o'])
    print('\nDONE: Analysis finished succesfully.\n')
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
    required_args = {'-f', '-r', '-d', '-m', '-o','-s'}
    arg_value_types = {'-f': str, '-r': str, '-d': str, '-m': str, '-o': str, '-D': int, '-q': int, '-c': float, '-i': float, '-g': str,'--adapter_sequence':str, '--adapter_sequence_2':str,'-s':str}
    min_arg_values = {'-D': 1, '-q': 0, '-c': 0, '-i': 0}
    max_arg_values = {'-c': 100, '-i': 100}
    default_arg_values = {'-D': 20, '-q': 30, '-c': 1.5, '-i': 90, '-g': 'yes','--adapter_sequence':'', '--adapter_sequence_2':''}
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

def cut_primers(output, consensus_seqs):
    '''Remove primers from consensus sequences before genotyping.'''
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

    consensus_name = list(seqs.keys())
    sequences = list(seqs.values())
    #consensus_seqs_no_primer = os.path.join(output, output + '_consensus_seqs_noprimers.fa')
    with open(consensus_seqs, 'w') as output_file:
        for index in range(len(consensus_name)):
            header = consensus_name[index]
            if '|core|' in header:
                sequences[index] = sequences[index][22:382]
            if '|ns5b|' in header:
                sequences[index] = sequences[index][22:364]
            output_file.write('>'+header + '\n')
            output_file.write(sequences[index]  + '\n')
    return consensus_seqs


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
    #terminal_command = f'samtools view -f 3 -F 2828 -q 30 -h {sam_out} | samtools sort -o {bam_out}'
    terminal_command = f'samtools view -f 1 -F 2828 -q 30 -h {sam_out} | samtools sort -o {bam_out}'
    error_msg = f'samtools view/sort terminated with errors while filtering and sorting mapped reads. Please refer to /{output}/logs/ for output logs.'
    stdout_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stdout.txt')
    stderr_file = os.path.join(output, 'logs', output + '_samtools_view_sort_stderr.txt')
    run(terminal_command, error_msg, stdout_file, stderr_file)

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
    blast_results = blast_results[blast_results['pident'] >= 90]
    #blast_results['end'] = blast_results.apply(lambda row: row['sseqid'].split('|')[1], axis=1)
    # Discard alignments below minimum identity threshold

    best_bitscores = blast_results[['qseqid','bitscore']].groupby(['qseqid']).head(10)
    blast_results = pd.merge(blast_results, best_bitscores, on=['qseqid', 'bitscore'])
    blast_results = blast_results.drop_duplicates()
    #blast_results1 = blast_results.sort_values(by=['bitscore'],ascending=False)
    blast_results.to_csv(os.path.join(output, output + '_genotype_calls.csv'), index=False)

    return blast_results

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
    #spades_out = os.path.join(output, output + '_spades_results')
    #sh.rmtree(spades_out)
    #blast_out = os.path.join(output, output + '_contigs_blast_results.tsv')
    #os.remove(blast_out)
    #contigs = os.path.join(output, output + '_contigs.fa')
    #os.remove(contigs)
    sam_out = os.path.join(output, output + '_alignment.sam')
    os.remove(sam_out)
    #ref_seqs = os.path.join(output, output + '_ref_seqs_for_mapping.fa')
    vcf_out = os.path.join(output, output + '_variants.vcf.gz')
    low_cov = os.path.join(output, output + '_low_cov.bed')
    #files = [ref_seqs + suffix for suffix in ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']]
    files = [vcf_out + suffix for suffix in ['', '.csi']]
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
