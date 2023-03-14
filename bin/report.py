#!/usr/bin/env python3

import pandas as pd
import re
import numpy as np
import argparse


def main(args):
    #consensus_report = pd.read_csv('combined_consensus_seqs_report.tsv', sep='\t')
    #genotype_calls = pd.read_csv('combined_genotype_calls.csv')
    #demix_report = pd.read_csv('combined_demix_report.csv')
    #qc_report = pd.read_csv('combined_qc_stats.csv')
    #basic_qc = pd.read_csv('../basic_qc/basic_qc_stats.csv')
    #abundance = pd.read_csv('../abundance_top_n/top_5_abundances_species.csv')
    consensus_report = pd.read_csv(args.consensus_report, sep='\t')
    genotype_calls = pd.read_csv(args.genotype_calls)
    demix_report = pd.read_csv(args.demix_report)
    qc_report = pd.read_csv(args.qc_report)
    basic_qc = pd.read_csv(args.basic_qc)
    abundance = pd.read_csv(args.abundance_top_n)

    with open(args.fastqlist) as f:
        fqlist = f.readlines()

    fqlist = [q.strip() for q in fqlist]

    fqlist_df = pd.DataFrame(fqlist, columns = ['sample_id'])

    #process consensus result
    consensus_report['sample_id'] = consensus_report['consensus_seq'].apply(lambda x: re.sub("^_R_","", x))
    consensus_report['sample_id'] = consensus_report['sample_id'].apply(lambda x: x.split('_')[0])

    consensus_report_reduced = consensus_report[["sample_id","amplicon","subtype","sequenced_bases"]]
    consensus_report_reduced =  consensus_report_reduced.astype(str)
    consensus_report_concat = consensus_report_reduced.groupby(['sample_id', 'amplicon']).agg({'subtype': [('subtype', '|'.join)],'sequenced_bases':[('sequenced_bases','|'.join)]}).reset_index()
    consensus_report_concat = pd.DataFrame(consensus_report_concat.values, columns=['sample_id','amplicon','subtype','sequenced_bases'])
    consensus_report_spread = consensus_report_concat.pivot(index=['sample_id'],columns='amplicon',values=['subtype','sequenced_bases']).reset_index()
    consensus_report_spread.columns = consensus_report_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    consensus_report_spread = consensus_report_spread.rename(columns={'_sample_id': 'sample_id'})
    #process qc results


    qc_report = pd.merge(consensus_report['consensus_seq'],qc_report, left_on="consensus_seq", right_on="contig",how='left')
    qc_report_reduced = qc_report[["sample_id","amplicon","total_mapped_bases","mean_coverage","std_coverage","proportion_genome_covered_over_20x"]]
    qc_report_reduced[['mean_coverage']] = round(qc_report_reduced[['mean_coverage']]).astype(int)
    qc_report_reduced[['std_coverage']] = round(qc_report_reduced[['std_coverage']]).astype(int)
    qc_report_reduced[['proportion_genome_covered_over_20x']] = round(qc_report_reduced[['proportion_genome_covered_over_20x']],2)
    qc_report_reduced =  qc_report_reduced.astype(str)
    qc_report_concat = qc_report_reduced.groupby(['sample_id', 'amplicon']).agg({'total_mapped_bases': [('total_mapped_bases', '|'.join)],'mean_coverage':[('mean_coverage','|'.join)],'std_coverage':[('std_coverage','|'.join)],'proportion_genome_covered_over_20x':[('proportion_genome_covered_over_20x','|'.join)]}).reset_index()
    qc_report_concat = pd.DataFrame(qc_report_concat.values, columns=['sample_id','amplicon','total_mapped_bases','mean_coverage',"std_coverage","proportion_genome_covered_over_20x"])
    qc_report_spread = qc_report_concat.pivot(index=['sample_id'],columns='amplicon',values=['total_mapped_bases','mean_coverage','std_coverage','proportion_genome_covered_over_20x']).reset_index()
    qc_report_spread.columns = qc_report_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    qc_report_spread = qc_report_spread.rename(columns={'_sample_id': 'sample_id'})

    #process genotype results
    genotype_calls['sample_id'] = genotype_calls['qseqid'].apply(lambda x: re.sub("^_R_","", x))
    genotype_calls['sample_id'] = genotype_calls['sample_id'].apply(lambda x: x.split('_')[0])
    genotype_calls['amplicon'] = genotype_calls['qseqid'].apply(lambda x: x.split('|')[1])
    genotype_calls['qseqid'] = genotype_calls['qseqid'].apply(lambda x: x.split('|')[0])

    genotype_calls_reduced = genotype_calls[['sample_id','qseqid', 'sseqid','amplicon','bitscore']]

    #top1 genotype
    best_bitscores = genotype_calls_reduced.groupby(['qseqid','amplicon']).max('bitcore').reset_index()
    best_gt = pd.merge(genotype_calls_reduced, best_bitscores, on=['qseqid', 'amplicon','bitscore'])

    #best_gt = best_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',np.mean)]}).reset_index()
    best_gt = best_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',list)]}).reset_index()
    best_gt = pd.DataFrame(best_gt.values, columns = ['sample_id','amplicon','top1_genotype','top1_bitscore'])
    best_gt_spread = best_gt.pivot(index=['sample_id'],columns='amplicon',values=['top1_genotype','top1_bitscore']).reset_index()
    best_gt_spread.columns = best_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    best_gt_spread = best_gt_spread.rename(columns={'_sample_id': 'sample_id'})

    #top2 genotype

    res = pd.merge(genotype_calls_reduced,best_bitscores,  on=['qseqid', 'amplicon','bitscore'], how='outer',indicator=True).query('_merge != "both"').drop(columns='_merge')
    second_bitscores = res.groupby(['qseqid','amplicon']).max('bitcore').reset_index()
    second_gt = pd.merge(res, second_bitscores, on=['qseqid', 'amplicon','bitscore'])

    second_gt = second_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',list)]}).reset_index()
    second_gt = pd.DataFrame(second_gt.values, columns = ['sample_id','amplicon','top2_genotype','top2_bitscore'])
    second_gt_spread = second_gt.pivot(index=['sample_id'],columns='amplicon',values=['top2_genotype','top2_bitscore']).reset_index()
    second_gt_spread.columns = second_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    second_gt_spread = second_gt_spread.rename(columns={'_sample_id': 'sample_id'})
    #merge 1 with consensus report
    merge0 = pd.merge(fqlist_df,basic_qc[['sample_id','total_bases']],on = 'sample_id', how='left')
    merge01 = pd.merge(merge0,abundance,on = 'sample_id', how='left')
    merge1 = pd.merge(merge01,consensus_report_spread, on = 'sample_id', how='left')
    merge2 = pd.merge(merge1, qc_report_spread, on='sample_id',how='left')
    merge3 = pd.merge(merge2, best_gt_spread, on='sample_id',how='left')
    merge4 = pd.merge(merge3, second_gt_spread, on='sample_id',how='left')
    merge5 = pd.merge(merge4, demix_report, on='sample_id',how='left')

    conditions = [
      (merge5['abundance_1_name'] == "Hepacivirus C"),
      (merge5['abundance_2_name'] == "Hepacivirus C"),
      (merge5['abundance_3_name'] == "Hepacivirus C"),
      (merge5['abundance_4_name'] == "Hepacivirus C"),
      (merge5['abundance_5_name'] == "Hepacivirus C")
    ]
    choices = [merge5['abundance_1_fraction_total_reads']* merge5['total_bases'],merge5['abundance_2_fraction_total_reads']*merge5['total_bases'],merge5['abundance_3_fraction_total_reads']* merge5['total_bases'],merge5['abundance_4_fraction_total_reads']* merge5['total_bases'],merge5['abundance_5_fraction_total_reads']* merge5['total_bases']]
    merge5['hcv_bases'] =  np.select(conditions,choices,default=0) 
    #print(merge5)
    #FLAGGING samples
    conditions = [
        #~(merge5['most_abundant_species_name'] == "Hepacivirus C"),#check if hcv is not the most abundant, even for negative controls
        (merge5['hcv_bases'] < 22000),
        (merge5['core_subtype'].isna() | merge5['ns5b_subtype'].isna()),
        ~(merge5['core_subtype'] == merge5['ns5b_subtype']) & (merge5['core_subtype'].notna() & merge5['ns5b_subtype'].notna()),
        (merge5['ns5b_sequenced_bases'].notna()) & (merge5['ns5b_sequenced_bases'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 300),
        (merge5['core_sequenced_bases'].notna()) & (merge5['core_sequenced_bases'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 300),
        (merge5['core_mean_coverage'].notna()) & (merge5['core_mean_coverage'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 20),
        (merge5['ns5b_mean_coverage'].notna()) & (merge5['ns5b_mean_coverage'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 20)
        
    ]

    choices = ['Check - low HCV content','Check - missing core/ns5b subtype', 'Check - mismatch core/ns5b subtypes','Check - ns5b sequenced bases','Check - core sequenced bases', 'Check - core mean coverage < 20', 'Check - ns5b mean coverage < 20']

    merge5['check'] = np.select(conditions,choices,default="PASS") 

    merge5.to_csv('run_summary_report.csv',index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--consensus_report')
    parser.add_argument('--genotype_calls')
    parser.add_argument('--demix_report')
    parser.add_argument('--qc_report')
    parser.add_argument('--basic_qc')
    parser.add_argument('--abundance_top_n')
    parser.add_argument('--fastqlist')
#    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
