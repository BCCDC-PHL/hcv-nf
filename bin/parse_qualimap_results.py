#!/usr/bin/env python

import argparse
import csv
import json
import sys
import pandas as pd

def main(args):

    with open(args.qualimap_bamqc_genome_results, 'r') as f:
    
        contigs =[]
     
        output_data = {}
        
        #output_data['genotype'] = genotype_name[0].split("|")[2]
        #output_data['segment'] = genotype_name[0].split("|")[1]
        for line in f:
            line = line.strip()
            if line.startswith('bam file'):
                sample_id = line.split('=')[1].strip().split('_')[0]
                output_data['sample_id'] = sample_id
            if line.startswith('median insert size'):
                median_insert_size = line.split('=')[1].strip().replace(',', '')
                output_data['median_insert_size'] = int(median_insert_size)
            #if line.startswith('mean coverageData'):
            #    mean_coverage = line.split('=')[1].strip().strip('X')
            #    output_data['mean_coverage'] = round(float(mean_coverage.replace(',','')), 2)
            #if line.startswith('std coverageData'):
            #    stdev_coverage = line.split('=')[1].strip().strip('X')
            #    output_data['stdev_coverage'] = round(float(stdev_coverage.replace(',','')), 2)
            if 'reference with a coverageData >= 5X' in line:
                proportion_genome_covered_over_5x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_5x'] = round(proportion_genome_covered_over_5x, 4)
            if 'reference with a coverageData >= 10X' in line:
                proportion_genome_covered_over_10x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_10x'] = round(proportion_genome_covered_over_10x, 4)
            if 'reference with a coverageData >= 20X' in line:
                proportion_genome_covered_over_20x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_20x'] = round(proportion_genome_covered_over_20x, 4)
            if 'reference with a coverageData >= 30X' in line:
                proportion_genome_covered_over_30x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_30x'] = round(proportion_genome_covered_over_30x, 4)
            if 'reference with a coverageData >= 40X' in line:
                proportion_genome_covered_over_40x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_40x'] = round(proportion_genome_covered_over_40x, 4)
            if 'reference with a coverageData >= 50X' in line:
                proportion_genome_covered_over_50x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_50x'] = round(proportion_genome_covered_over_50x, 4)
            if '|' in line:
                contigs.append(line)


    data = []
    for i in contigs:
        contig = i.strip()
        
        if len(contig) > 0:
            data.append([contig.split('\t')[0].split('|')[0].split('_')[0],contig.split('\t')[0].split('|')[0], contig.split('\t')[0].split('|')[1],contig.split('\t')[0].split('|')[2],contig.split('\t')[1],contig.split('\t')[2],contig.split('\t')[3],contig.split('\t')[4]])
            
    data_df = pd.DataFrame(data)
    data_df.columns = ['sample_id', 'contig','segment','subtype','seq_length','total_mapped_bases','mean_coverage','std_coverage']
    output_data_df = pd.DataFrame(output_data, index=[0])

    outdata = data_df.merge(output_data_df,how="left",on='sample_id')
    #output_data = pd.DataFrame(output_data, index=[0])
    outdata.to_csv(args.outfile,index=False)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qualimap_bamqc_genome_results')
    parser.add_argument('-o','--outfile')
#    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)

