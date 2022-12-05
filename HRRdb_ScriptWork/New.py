#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:37:53 2022

@author: nilesh@4basecare.com
"""

#### Package Import
import pandas as pd
import os
from pathlib import Path
#### Package Import

#%%  Edit Variables Here

dirpath = '/home/bioinfo/git/HRRdb/HRRdb_ScriptWork/'
hrr_genes = pd.read_csv("/home/bioinfo/git/HRRdb/HRR_genes.txt", sep='/t')
chunk_size = 3

#%% Listing of Files and Chunk Creation

file_list = os.listdir(dirpath)
file_list.remove('HRRdb_filemaker.py')
file_list.remove('New.py')

chunked_list = [file_list[i:i+chunk_size] for i in range(0, len(file_list), chunk_size)]
#print(chunked_list)

#%% Copying of Files to Folders created based on Chunked_List Size

dt=dict(enumerate(chunked_list))

for key in dt:
  print(key)
  os.mkdir(str(key))
  
dt_list = []
for key, values in dt.items():
    for value in values:
        dt_list.append([key, value])

dt_list = [[k,v] for k, values in dt.items() for v in values]

for a,b in dt_list:
    os.system('cp '+dirpath+str(b)+' '+str(a)+'/'+str(b))
    

#%% Working on Files in each Chunk Folder

for key in dt:
    os.chdir(str(key))
    #keydir_list = os.listdir(os.chdir(str(key)))
    key_files = sorted(Path('./').glob('*.csv'))
    
    samples = pd.DataFrame()
    for file in key_files:
        #f_path = dirpath + files
        file_df = pd.read_csv(file, sep=',')
        file_df['Sample_ID'] = file.name.strip('.csv')
        samples = samples.append(file_df)
        column_to_move = samples.pop('Sample_ID')
        samples.insert(0, "Sample_ID", column_to_move)
        hrr_df = pd.merge(samples, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
        del hrr_df["HRR_Genes"]
        hrr_df.to_csv('./'+str(key)+'_hrr_df.csv', index = False)
    os.chdir('../')
