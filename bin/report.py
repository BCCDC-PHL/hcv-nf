#!/usr/bin/env python3

import pandas as pd
import re
import numpy as np
import argparse


def main(args):

    consensus_report = pd.read_csv(args.consensus_report, sep='\t')
    demix_report = pd.read_csv(args.demix_report)
    qc_report = pd.read_csv(args.qc_report)
    mapped_reads_counts = pd.read_csv(args.reads_counts)
    genotype_calls = pd.read_csv(args.genotype_nt)
    #genotype_calls = pd.read_csv("combined_genotype_calls.csv")

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
    genotype_calls['sample_id'] = genotype_calls['query_seq_id'].apply(lambda x: re.sub("^_R_","", x))
    genotype_calls['sample_id'] = genotype_calls['sample_id'].apply(lambda x: x.split('_')[0])
    genotype_calls['amplicon'] = genotype_calls['query_seq_id'].apply(lambda x: x.split('|')[1])
    genotype_calls['query_seq_id'] = genotype_calls['query_seq_id'].apply(lambda x: x.split('|')[0])

    gt_counts = genotype_calls.groupby(['sample_id','subject_names']).size().reset_index(name = 'counts')
    gt_counts = gt_counts[gt_counts['subject_names'].str.contains('genotype') | gt_counts['subject_names'].str.contains('subtype') ]

    gt_counts['nt_genotypes'] = gt_counts['subject_names'] + '(n=' + gt_counts['counts'].astype(str) + ')'
    gt_counts['nt_genotypes'] = gt_counts['nt_genotypes'].apply(lambda x: re.sub(r'^.*?genotype', '', x).strip()).apply(lambda x: re.sub(r'^.*?subtype', '', x).strip() )
    gt_nt = gt_counts.groupby("sample_id")["nt_genotypes"].apply(" + ".join).reset_index(name = "nt_genotypes")


    #genotype_calls_reduced = genotype_calls[['sample_id','qseqid', 'sseqid','amplicon','bitscore']]

    #top1 genotype
    #best_bitscores = genotype_calls_reduced.groupby(['qseqid','amplicon']).max('bitcore').reset_index()
    #best_gt = pd.merge(genotype_calls_reduced, best_bitscores, on=['qseqid', 'amplicon','bitscore'])

    #best_gt = best_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',list)]}).reset_index()
    #best_gt = pd.DataFrame(best_gt.values, columns = ['sample_id','amplicon','top1_genotype','top1_bitscore'])
    #best_gt_spread = best_gt.pivot(index=['sample_id'],columns='amplicon',values=['top1_genotype','top1_bitscore']).reset_index()
    #best_gt_spread.columns = best_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    #best_gt_spread = best_gt_spread.rename(columns={'_sample_id': 'sample_id'})

    #top2 genotype

    #res = pd.merge(genotype_calls_reduced,best_bitscores,  on=['qseqid', 'amplicon','bitscore'], how='outer',indicator=True).query('_merge != "both"').drop(columns='_merge')
    #second_bitscores = res.groupby(['qseqid','amplicon']).max('bitcore').reset_index()
    #second_gt = pd.merge(res, second_bitscores, on=['qseqid', 'amplicon','bitscore'])

    #second_gt = second_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',list)]}).reset_index()
    #second_gt = pd.DataFrame(second_gt.values, columns = ['sample_id','amplicon','top2_genotype','top2_bitscore'])
    #second_gt_spread = second_gt.pivot(index=['sample_id'],columns='amplicon',values=['top2_genotype','top2_bitscore']).reset_index()
    #second_gt_spread.columns = second_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
    #second_gt_spread = second_gt_spread.rename(columns={'_sample_id': 'sample_id'})
    #merge 1 with consensus report
    #merge0 = pd.merge(fqlist_df,basic_qc[['sample_id','total_bases']],on = 'sample_id', how='left')
    #merge01 = pd.merge(merge0,abundance,on = 'sample_id', how='left')
    merge1 = pd.merge(fqlist_df,consensus_report_spread, on = 'sample_id', how='left')
    merge2 = pd.merge(merge1, qc_report_spread, on='sample_id',how='left')
  
    merge5 = pd.merge(merge2, demix_report, on='sample_id',how='left')
    merge6 = pd.merge(merge5, mapped_reads_counts, on='sample_id',how='left')
    merge7 = pd.merge(merge6, gt_nt, on='sample_id', how = 'left')

    conditions = [

        ((merge7['core_subtype'].isna()) & (merge7['core_mapped_reads'] < 50)),
        ((merge7['ns5b_subtype'].isna()) & (merge7['ns5b_mapped_reads'] < 50)),
        (merge7['core_subtype'].isna() | merge7['ns5b_subtype'].isna()),
        (merge7['ns5b_sequenced_bases'].notna()) & (merge7['ns5b_sequenced_bases'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 300),
        (merge7['core_sequenced_bases'].notna()) & (merge7['core_sequenced_bases'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 300),
        (merge7['core_mean_coverage'].notna()) & (merge7['core_mean_coverage'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 20),
        (merge7['ns5b_mean_coverage'].notna()) & (merge7['ns5b_mean_coverage'].astype(str).apply(lambda x: min(x.split('|'))).replace('nan',np.nan).fillna(0).astype(int) < 35),
        ~(merge7['core_subtype'] == merge7['ns5b_subtype']) & (merge7['core_subtype'].notna() & merge7['ns5b_subtype'].notna())
        
    ]

    choices = ['Check - core low mapped reads','Check - ns5b low mapped reads', 'Check - missing core/ns5b subtype', 'Check - ns5b sequenced bases','Check - core sequenced bases', 'Check - core mean coverage < 20', 'Check - ns5b mean coverage < 35','Check - mismatch core/ns5b subtypes']

    merge7['check'] = np.select(conditions,choices,default="PASS") 

    merge7.to_csv(args.prefix + '_run_summary_report.csv',index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--consensus_report')
    parser.add_argument('--demix_report')
    parser.add_argument('--qc_report')
    parser.add_argument('--reads_counts')
    parser.add_argument('--fastqlist')
    parser.add_argument('--genotype_nt')
    parser.add_argument('--prefix')
    args = parser.parse_args()
    main(args)
