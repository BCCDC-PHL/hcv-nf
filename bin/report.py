#!/usr/bin/env python3

import pandas as pd
import re
import numpy as np

consensus_report = pd.read_csv('combined_consensus_seqs_report.tsv', sep='\t')
genotype_calls = pd.read_csv('combined_genotype_calls.csv')
demix_report = pd.read_csv('combined_demix_report.csv')
qc_report = pd.read_csv('combined_qc_stats.csv')

with open('fastqlist') as f:
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


qc_report_reduced = qc_report[["sample_id","amplicon","total_mapped_bases","mean_coverage","std_coverage","proportion_genome_covered_over_20x"]]
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

genotype_calls_reduced = genotype_calls[['sample_id','sseqid','amplicon','bitscore']]

#top1 genotype
best_bitscores = genotype_calls_reduced.groupby(['sample_id','amplicon']).max('bitcore').reset_index()
best_gt = pd.merge(genotype_calls_reduced, best_bitscores, on=['sample_id', 'amplicon','bitscore'])

best_gt = best_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',np.mean)]}).reset_index()
best_gt = pd.DataFrame(best_gt.values, columns = ['sample_id','amplicon','top1_genotype','top1_bitscore'])
best_gt_spread = best_gt.pivot(index=['sample_id'],columns='amplicon',values=['top1_genotype','top1_bitscore']).reset_index()
best_gt_spread.columns = best_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
best_gt_spread = best_gt_spread.rename(columns={'_sample_id': 'sample_id'})

#top2 genotype

res = pd.merge(genotype_calls_reduced,best_bitscores,  on=['sample_id', 'amplicon','bitscore'], how='outer',indicator=True).query('_merge != "both"').drop(columns='_merge')
second_bitscores = res.groupby(['sample_id','amplicon']).max('bitcore').reset_index()
second_gt = pd.merge(res, second_bitscores, on=['sample_id', 'amplicon','bitscore'])

second_gt = second_gt.groupby(['sample_id', 'amplicon']).agg({'sseqid': [('sseqid', ','.join)],'bitscore':[('bitscore',np.mean)]}).reset_index()
second_gt = pd.DataFrame(second_gt.values, columns = ['sample_id','amplicon','top2_genotype','top2_bitscore'])
second_gt_spread = second_gt.pivot(index=['sample_id'],columns='amplicon',values=['top2_genotype','top2_bitscore']).reset_index()
second_gt_spread.columns = second_gt_spread.columns.map(lambda index: f'{index[1]}_{index[0]}')
second_gt_spread = second_gt_spread.rename(columns={'_sample_id': 'sample_id'})
#merge 1 with consensus report
merge1 = pd.merge(fqlist_df,consensus_report_spread, on = 'sample_id', how='left')
merge2 = pd.merge(merge1, qc_report_spread, on='sample_id',how='left')
merge3 = pd.merge(merge2, best_gt_spread, on='sample_id',how='left')
merge4 = pd.merge(merge3, second_gt_spread, on='sample_id',how='left')
merge5 = pd.merge(merge4, demix_report, on='sample_id',how='left')

merge5.to_csv('run_summary_report.csv',index=False)