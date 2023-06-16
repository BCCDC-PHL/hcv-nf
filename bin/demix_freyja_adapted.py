#!/usr/bin/env python3

import pandas as pd
import argparse
import os
from freyja.sample_deconv import build_mix_and_depth_arrays,solve_demixing_problem,map_to_constellation


def reindex_dfs(df_barcodes, mix, depths):
    # reindex everything to match across the dfs
    df_barcodes = df_barcodes.reindex(sorted(df_barcodes.columns), axis=1)
    mix = mix.reindex(df_barcodes.columns).fillna(0.)

    mix_as_set = set(mix.index)
    # dropping extra sequencing depth info we don't need
    depths = depths.drop(labels=[m0 for m0 in df_barcodes.columns
                                 if m0 not in mix_as_set])
    depths = depths.reindex(df_barcodes.columns).fillna(0.)
    return df_barcodes, mix, depths

def demix(variants, depths, output, eps, barcodes, covcut):
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '-1':
        df_barcodes = pd.read_csv(barcodes, index_col=0)
    else:
        df_barcodes = pd.read_csv(os.path.join(locDir,'data/barcodes_core.csv'), index_col=0)

    muts = list(df_barcodes.columns)
    #mapDict = buildLineageMap(meta)
    #mapDict = buildLineageMap(meta)
    #create mapdict
    dicts = {}
    lineages = df_barcodes.index
    for i in lineages:
        dicts[i] = i.split('_')[0]

    print('building mix/depth matrices')
    # assemble data from (possibly) mixed samples
    mix, depths_, cov = build_mix_and_depth_arrays(variants, depths, muts,
                                                   covcut)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                               mix,
                                                               depths_,
                                                               eps)
    # merge intra-lineage diversity if multiple hits.
    if len(set(sample_strains)) < len(sample_strains):
        localDict = {}
        for jj, lin in enumerate(sample_strains):
            if lin not in localDict.keys():
                localDict[lin] = abundances[jj]
            else:
                localDict[lin] += abundances[jj]
        # ensure descending order
        localDict = dict(sorted(localDict.items(),
                                key=lambda x: x[1],
                                reverse=True))
        sample_strains = list(localDict.keys())
        abundances = list(localDict.values())
    
    localDict = map_to_constellation(sample_strains, abundances, dicts)
    
    # assemble into series and write.
    sols_df = pd.Series(data=(localDict,sample_strains, abundances,
                              error, cov),
                        index=['summarized','lineages',
                        'abundances', 'resid', 'coverage'],
                        name=mix.name)
    # convert lineage/abundance readouts to single line strings
    sols_df['lineages'] = ' '.join(sols_df['lineages'])
    sols_df['abundances'] = ['%.8f' % ab for ab in sols_df['abundances']]
    
    sols_df.to_csv(output, sep='\t')
    
    reformated = []
    for i in localDict:
        i_ls = list(i)
        i_ls[1] = '%.4f' % i_ls[1]
        i = tuple(i_ls)
        reformated.append(i)


    result = pd.DataFrame(data=[[args.sample,reformated]], columns=['sample_id','demix_results_core'])
    result.to_csv(args.sample+'_demix.csv',index = False)

def main(args):
    demix(args.variants, args.depths, args.output, args.eps, args.barcodes, args.covcut)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('variants')
    parser.add_argument('depths')
    parser.add_argument('--eps', default=1e-3, help='minimum abundance to include')
    parser.add_argument('--barcodes', default='-1', help='custom barcode file')
    parser.add_argument('--output', default='demixing_result.csv', help='Output file')
    parser.add_argument('--sample', default='', help='Output file')
    parser.add_argument('--covcut', default=10, help='depth cutoff for\
                                            coverage estimate')
    args = parser.parse_args()
    main(args)
