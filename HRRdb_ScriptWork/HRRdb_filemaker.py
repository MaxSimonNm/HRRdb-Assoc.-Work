#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:36:32 2022

@author: nilesh@4basecare.com
"""
#### Package Import
"""
Created on Tue Nov 15 17:36:32 2022
@author: nilesh@4basecare.com
"""
#### Package Import
import pandas as pd
import numpy as np
import os
#### Package Import

wd = os.getcwd()

files = os.listdir(wd)
files.remove('HRRdb_filemaker.py')


filegrp=[]
it = iter(files)
for x, y in zip(it, it):
    print (x, y)
    file_tuple = (x,y)
    filegrp.append(file_tuple)
    
hrr_genes = pd.read_csv("/home/max/Github/HRRdb/HRR_genes.txt", sep='/t')

#hrr_df = pd.merge(df_all, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')

#hrr_df.pop("HRR_Genes")
#del hrr_df["HRR_Genes"]  #Deleting the HRR_Genes col

for (j, k) in filegrp:
    #print(j)
    #print(k)
    df1 = pd.read_csv(j)
    #f = j+".csv"
    #df1.to_csv(f, index=False)
    
    hrr_df1 = pd.merge(df1, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
    #hrr_df.pop("HRR_Genes")
    del hrr_df1["HRR_Genes"]  #Deleting the HRR_Genes col
    
    hrr_df1["Sample_ID"] = np.nan
    hrr_df1["Sample_ID"] = hrr_df1["Sample_ID"].fillna(j)
    
    column_to_move = hrr_df1.pop("Sample_ID") # insert column with insert(location, column_name, column_value)
    hrr_df1.insert(0, "Sample_ID", column_to_move)
        
    
    df2 = pd.read_csv(k)
    #g = k+".csv"
    #df2.to_csv(g, index=False)
    
    hrr_df2 = pd.merge(df2, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
    #hrr_df.pop("HRR_Genes")
    del hrr_df2["HRR_Genes"]  #Deleting the HRR_Genes col
    
    hrr_df2["Sample_ID"] = np.nan
    hrr_df2["Sample_ID"] = hrr_df2["Sample_ID"].fillna(k)
    
    column_to_move = hrr_df2.pop("Sample_ID") # insert column with insert(location, column_name, column_value)
    hrr_df2.insert(0, "Sample_ID", column_to_move)
    
    #~ Export
    
    filename = j.replace('_output.csv', '') + "_"+ k.replace('_output.csv', '') + ".csv"
        
    frames = [hrr_df1, hrr_df2]
    result = pd.concat(frames)
    result.to_csv(filename, index=False)
    
      











'''
#%% Module 1 Merging Files
source_files = sorted(Path('/home/max/Github/HRRdb/HRRdb_ScriptWork/').glob('*.csv'))
#%%
dataframes =[]
for group in chunker(source_files, 2):
    print(group)
    for file in group:
        print(file)
        df = pd.read_csv(file, sep=',') # additional arguments up to your needs
        df['Sample_ID'] = file.name
        column_to_move = df.pop('Sample_ID') # insert column with insert(location, column_name, column_value)
        df.insert(0, "Sample_ID", column_to_move)
        
        hrr_genes = pd.read_csv("/home/max/Github/HRRdb/HRR_genes.txt", sep='/t')
        hrr_df = pd.merge(df, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
        del hrr_df["HRR_Genes"]  #Deleting the HRR_Genes col
        hrr_df.to_csv(nameoffile)
        #[(file1, file2), (file3,file4)]
        
        
        
        
        #dataframes.append(df)
        
      
        
#%%
df_all = pd.concat(dataframes)
column_to_move = df_all.pop('Sample_ID') # insert column with insert(location, column_name, column_value)
df_all.insert(0, "Sample_ID", column_to_move)
#~ Export
#df_all.to_csv('/home/bioinfo/Nilesh/FE_merged/merged.csv', index = False)
#%%Module 2 Filtering Files based on HRR Genes
hrr_genes = pd.read_csv("/home/max/Github/HRRdb/HRR_genes.txt", sep='/t')
hrr_df = pd.merge(df_all, hrr_genes, how="inner", left_on='Ref.Gene', right_on='HRR_Genes')
#hrr_df.pop("HRR_Genes")
del hrr_df["HRR_Genes"]  #Deleting the HRR_Genes col
#~ Export
#hrr_df.to_csv('/home/max/Github/HRRdb/HRRdb_ScriptWork/df_new.csv', index = False)
'''