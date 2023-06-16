#!/usr/bin/env python3
import argparse
import os
import subprocess as sp
import shutil as sh


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

def main(args):

    #read in best contigs
    output = args.sample
    with open(args.input, 'r') as input_file:
    #with open("ns5b_contig.fa", 'r') as input_file:
        contig_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                contig_seqs[header] = ''
            else:
                contig_seqs[header] += line.strip()

    names = list(contig_seqs)

    ref_names = [x.split('|')[3] + '_' + x.split('|')[2] for x in names]
    #read in best ref
    with open(args.ref, 'r') as input_file:
    #with open("best_ns5b_ref.fa", 'r') as input_file:  
        ref_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>').split(' ')[0]
                ref_seqs[header] = ''
            else:
                ref_seqs[header] += line.strip()


    ref_seqs_for_mapping = []
    for n in range(len(ref_names)):      
        filename = names[n].split('|')[0] +  '_' + names[n].split('|')[1] + '_' + ref_names[n]  + '_contig_ref.fa'
        prefix = names[n].split('|')[0] +  '_' + names[n].split('|')[1] + '_' + ref_names[n] 
        with open(filename,'w') as output_file:
            output_file.write('>'+ ref_names[n] + '\n')
            output_file.write(ref_seqs[ref_names[n]] + '\n')
            output_file.write('>'+ names[n] + '\n')
            output_file.write(contig_seqs[names[n]] + '\n')
        terminal_command = (f'mafft --reorder --adjustdirection --anysymbol --thread 4 --auto {filename} > {prefix}_mafft.fa')      
        error_msg = f'mafft terminated with errors while aligning contigs to ref seqs. Please refer to /logs/ for output logs.'
        stdout_file = None
        stderr_file = os.path.join(output + '_mafft_stderr.txt')
        run(terminal_command, error_msg, stdout_file, stderr_file)
        with open(prefix+'_mafft.fa', 'r') as input_file:  
            ref_seq_for_mapping = {}
            for line in input_file:
                if line[0] == '>':
                    header = line.strip().lstrip('>').split(' ')[0]
                    ref_seq_for_mapping[header] = ''
                else:
                    ref_seq_for_mapping[header] += line.strip()     
        names_n = list(ref_seq_for_mapping)
        ref_seq = list(ref_seq_for_mapping[names_n[0]])
        contig_seq = list(ref_seq_for_mapping[names_n[1]])
        index_to_remove = []
        for s in range(len(ref_seq)):
            if (ref_seq[s] == '-') & ~(contig_seq[s] == '-'):
                contig_seq[s] = 'Z'
        contig_seq = ''.join(contig_seq).replace('Z','')
        contig_seq = contig_seq.replace('-','N')
        ref_seqs_for_mapping.append({names_n[1] : contig_seq })

 

    #ref_seq = list(contig_seqs[names[0]])

    #contig_seq = list(contig_seqs[names[1]])

    #index_to_remove = []
    #for s in range(len(ref_seq)):
    #    if (ref_seq[s] == '-') & ~(contig_seq[s] == '-'):
    #        contig_seq[s] = 'Z'

    #contig_seq = ''.join(contig_seq).replace('Z','')
    #contig_seq = contig_seq.replace('-','N')
    names_ref = list(ref_seqs_for_mapping)


    with open(args.outfile, 'w') as output_file:
        #contig_counter = 1
    #with open("ref_seqs_for_mapping_ns5b.fa", 'w') as output_file:
        for index in range(len(ref_seqs_for_mapping)):
            header = list(ref_seqs_for_mapping[index])[0]
            output_file.write('>'+header + '\n')
            output_file.write(ref_seqs_for_mapping[index][header] + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input')
    parser.add_argument('-o','--outfile')
    parser.add_argument('-s','--sample')
    parser.add_argument('-r','--ref')
#    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)

