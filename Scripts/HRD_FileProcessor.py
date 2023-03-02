#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:37:53 2022

@author: nilesh@4basecare.com
"""

#### Package Import
import warnings
import pandas as pd
import os
from pathlib import Path
import time
import telegram_send
####

warnings.filterwarnings("ignore")
start_time = time.time()

#%%  Edit Variables Here

#Folder containing MultiAnno Merged Output of Sample's VCF
dirpath = '/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/s4_Outside_1-17_Som_MultAnno/22_Outside_1-17_Som_MulAnno_merged/'

#Folder Path for Saving the Chunks
savepath = '/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/s4_Outside_1-17_Som_MultAnno/22_Outside_1-17_Som_HRD-Out/'

#Path for HRR Genes file
hrr_genes = pd.read_csv("/home/bioinfo/Nilesh/HRRdb_Samples/Scripts/HRR_genes.txt", sep='/t')

#Set Number of files per chunk.
##If total files is 5, and chunk size is 3,
##then 2 folders will be made,
##1st folder will have 3 files, 2nd one will have 2 files
chunk_size = 50
                
#%% Listing of Files and Chunk Creation based on Chunk Size

file_list = os.listdir(dirpath)

#Remove files that are not Samples
#file_list.remove('HRD_FileProcessor.py') 
#file_list.remove('save')

chunked_list = [file_list[i:i+chunk_size] for i in range(0, len(file_list), chunk_size)]
#print(chunked_list)

#%% Copying of Files to Folders created based on Chunked_List Size

dt=dict(enumerate(chunked_list))  #Nested List converted to Dictionary

for keys in dt:                   #Folder Creation based on number of Chunk keys
  print(keys)
  os.mkdir(savepath+str(keys))
  
dt_list = [[k,v] for k, values in dt.items() for v in values]
  
#dt_list = []                     #For Loop Expansion of above List Comprehension for Understanding
#for keys, values in dt.items():  #for making file list per chunk used later for copying
#   for value in values:
#       dt_list.append([keys, value])


for a,b in dt_list:               #Copying files to their specific chunk
    os.system('cp '+dirpath+str(b)+' '+savepath+str(a)+'/'+str(b))
    

#%% Working on Files in each Chunk Folder

for keys in dt:                   #Looping through each folder created based on keys
    os.chdir(savepath+str(keys))
    key_files = sorted(Path('./').glob('*.csv'))
    
    samples = pd.DataFrame()
    for file in key_files:
        file_df = pd.read_csv(file, sep=',')
        file_df['Sample_ID'] = file.name.strip('.csv')
        samples = samples.append(file_df)
        column_to_move = samples.pop('Sample_ID')
        samples.insert(0, "Sample_ID", column_to_move)
        #print(samples)


        
        # filter for rows containing "exonic" and "splicing", but not "ncR"
        samples_filtered = samples[samples["Func.ensGene"].str.contains("exonic|splicing|exonic;splicing") & ~samples["Func.ensGene"].str.contains("ncR")]
       
      
        # filter out rows containing the string "synonymous SNV"
        samples_filtered = samples_filtered[samples_filtered['ExonicFunc.ensGene'] != 'synonymous SNV']



        # filtering out pop freq columns for values less than 0.01
        popfreqs=['esp6500siv2_all','ExAC_ALL','ExAC_SAS','AF','AF_sas','1000g2015aug_all','1000g2015aug_SAS']
        
        for col in popfreqs:
            samples_filtered[col] = samples_filtered[col].replace('.',0).fillna(0) 
            samples_filtered[col] = samples_filtered[col].astype(float).round(4)
            samples_filtered=samples_filtered[samples_filtered[col]<= 0.01]
            
    
   
        #HRR Gene Filtration
        hrr_df = pd.merge(samples_filtered, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
        #print(hrr_df)
        del hrr_df["HRR_Genes"]
        hrr_df.to_csv('./'+str(keys)+'_hrr_df.csv', index = False)
    os.chdir('../')
#%%    
time_taken = ("Time taken: {:.2f} seconds".format(time.time() - start_time))
os.system("notify-send 'HRD File Processor' 'Process Finished'")    
msg = f"HRD File Processor finished @Nilesh_Mukherjee | {time_taken}"
print(msg) #for stdout and telegram-send --stdin
telegram_send.send(messages=[msg]) #for telegram notification
