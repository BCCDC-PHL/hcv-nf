#!/usr/bin/env python3

import sys as sys
import os as os
import subprocess as sp
import shutil as sh
import numpy as np
import pandas as pd


blast_results = pd.read_csv("A3-hcv220223CE-corens5b_filtered_blast_results.csv")
counts = blast_results['subtype'].value_counts().reset_index()

counts=counts.rename(columns = {'subtype' : 'counts','index':'subtype'})
total=counts[['counts']].sum()
counts['prop'] = counts['counts'].apply(lambda x: x/total)


mean_bitscore = blast_results[['subtype','bitscore']].groupby('subtype').mean().reset_index()

bit_df = pd.merge(counts,mean_bitscore,how ='left', on='subtype')
bit_df=bit_df.rename(columns = {'bitscore': 'mean_bitscore'})