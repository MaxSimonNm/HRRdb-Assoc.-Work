#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:36:32 2022

@author: nilesh@4basecare.com
"""
#### Package Import
import pandas as pd
from pathlib import Path
#import os as os
#### Package Import

####
source_files = sorted(Path('/home/bioinfo/Nilesh/FE_merged/').glob('*.csv'))

dataframes =[]
for file in source_files:
    df = pd.read_csv(file, sep=',') # additional arguments up to your needs
    df['Sample_ID'] = file.name
    dataframes.append(df)

df_all = pd.concat(dataframes)

column_to_move = df_all.pop('Sample_ID') # insert column with insert(location, column_name, column_value)

df_all.insert(0, "Sample_ID", column_to_move)

#df_all.to_csv('/home/bioinfo/Nilesh/FE_merged/merged.csv', index = False)

