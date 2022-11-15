#!/usr/bin/env python3
import argparse

def main(args):

    with open(args.input, 'r') as input_file:
        contig_seqs = {}
        for line in input_file:
            if line[0] == '>':
                header = line.strip().lstrip('>')
                contig_seqs[header] = ''
            else:
                contig_seqs[header] += line.strip()

    names = list(contig_seqs)

    ref_seq = list(contig_seqs[names[0]])

    contig_seq = list(contig_seqs[names[1]])

    index_to_remove = []
    for s in range(len(ref_seq)):
        if (ref_seq[s] == '-') & ~(contig_seq[s] == '-'):
            contig_seq[s] = 'Z'

    contig_seq = ''.join(contig_seq).replace('Z','')
    contig_seq = contig_seq.replace('-','N')



    with open(args.outfile, 'w') as output_file:
        #contig_counter = 1
        #for index in blast_results.index:
        header = '>'+names[1]
        output_file.write(header + '\n')
        output_file.write(contig_seq + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input')
    parser.add_argument('-o','--outfile')
#    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)

