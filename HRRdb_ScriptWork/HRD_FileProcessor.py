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
chunk_size = 3  #Set Number of files per chunk.
                #If total files is 5, and chunk size is 3,
                #then 2 folders will be made, 
                #1st folder will have 3 files, 2nd one will have 2 files

#%% Listing of Files and Chunk Creation based on Chunk Size

file_list = os.listdir(dirpath)

#Remove files that are not Samples
file_list.remove('HRRdb_filemaker.py') 
file_list.remove('New.py')

chunked_list = [file_list[i:i+chunk_size] for i in range(0, len(file_list), chunk_size)]
#print(chunked_list)

#%% Copying of Files to Folders created based on Chunked_List Size

dt=dict(enumerate(chunked_list))  #Nested List converted to Dictionary

for keys in dt:                   #Folder Creation based on number of Chunk keys
  print(keys)
  os.mkdir(str(keys))
  
dt_list = [[k,v] for k, values in dt.items() for v in values]
  
#dt_list = []                     #For Loop Expansion for above List Comprehension
#for keys, values in dt.items():  #for making file list per chunk used later for copying
#   for value in values:
#       dt_list.append([keys, value])


for a,b in dt_list:               #Copying files to their specific chunk
    os.system('cp '+dirpath+str(b)+' '+str(a)+'/'+str(b))
    

#%% Working on Files in each Chunk Folder

for keys in dt:                   #Looping through each folder created based on keys
    os.chdir(str(keys))
    key_files = sorted(Path('./').glob('*.csv'))
    
    samples = pd.DataFrame()
    for file in key_files:
        file_df = pd.read_csv(file, sep=',')
        file_df['Sample_ID'] = file.name.strip('.csv')
        samples = samples.append(file_df)
        column_to_move = samples.pop('Sample_ID')
        samples.insert(0, "Sample_ID", column_to_move)
        #print(samples)
        hrr_df = pd.merge(samples, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
        #print(hrr_df)
        del hrr_df["HRR_Genes"]
        hrr_df.to_csv('./'+str(keys)+'_hrr_df.csv', index = False)
    os.chdir('../')

os.system("notify-send 'HRD File Processor' 'Process Finished'")    
    
