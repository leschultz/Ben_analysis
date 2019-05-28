#!/usr/bin/env python3

'''
Construct APD dataframe of mean values.
'''

import pandas as pd

import sys
import os

df_file = sys.argv[1]  # File containing all Tg information
export_dir = sys.argv[2]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

df = pd.read_csv(df_file)
df = df.drop(['job'], axis=1)

cols = ['system', 'composition', 'hold_temperature']
groups = df.groupby(cols)

mean = groups.mean().add_suffix('_mean').reset_index()
sem = groups.sem().add_suffix('_sem').reset_index()
count = groups.count().add_suffix('_count').reset_index()

df = mean.merge(sem)
df = df.merge(count)

df.to_csv(os.path.join(export_dir, 'variance_mean_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'variance_mean_df.html'), index=False)
