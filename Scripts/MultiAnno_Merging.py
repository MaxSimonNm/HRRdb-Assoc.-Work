#!/usr/bin/env python3
import os
import pandas as pd
import warnings
import time
import telegram_send

start_time = time.time()

dirpath= '/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/s4_Outside_1-17_Som_MultAnno'
test='HRR'
sample_type='DNA [FE,FFPE]' #"DNA [Blood]"

#importing external files for filter engine
collist= pd.read_csv("/home/bioinfo/nishtha/TMB_17Aug/filter/columns39.csv")
canonical = pd.read_excel("/home/bioinfo/nishtha/TMB_17Aug/filter/canonical.xlsx", sheet_name=0, mangle_dupe_cols=True, engine='openpyxl')
genes= pd.read_csv("/home/bioinfo/nishtha/TMB_17Aug/filter/genelist.csv")
cohort4= pd.read_csv("/home/bioinfo/nishtha/TMB_17Aug/filter/4basecare-germline-cohort.csv")

testgenes= list(genes[test].dropna())
testgenes=[g.upper() for g in testgenes]

#making FE_merged folders in the destination dir
folders= os.listdir(dirpath)

if 'FE_merged' in folders:
    os.system('rm -r ' + dirpath + '/FE_merged')
folders= os.listdir(dirpath)    
os.system("mkdir " + dirpath + "/FE_merged")
warnings.filterwarnings("ignore")

        
tmb_filtered_df= pd.DataFrame(columns=['samplename','total_var','after pass','after allele freq', 'after mq', 'after dp', 'after exonic', 'after synony', 'after pop_freq', 'after allele freq2', 'after cohort4', 'tmb'])

#%%

#processing every sample in folder one by one
for f in folders:
    num=folders.index(f)
    f_path= dirpath + "/" + f
    print(f_path)
    files= os.listdir(f_path)
    cancer= [obj for obj in files if 'cancervar.hg19_multianno.txt.cancervar' in obj]
    annovar= [obj for obj in files if '_out.hg19_multianno' in obj]
    vcf= [obj for obj in files if 'final.tab' in obj]
    #print(f)

    #locations of different files
    cancerloc= dirpath + "/" + f + "/" + cancer[0]
    vcfloc= dirpath + "/" + f + "/" + vcf[0]
    annoloc= dirpath + "/" + f + "/" + annovar[0]

    #Detecting dragen 3.6 or 3.9 & reading columns file
    vcfdetect=pd.read_csv(vcfloc, sep="/t")
    if 'INFO:hotspot' in vcfdetect.columns:
        collist= pd.read_csv("/home/bioinfo/Nilesh/HRRdb_Samples/Scripts/columns39.csv")
    else:
        collist= pd.read_csv("/home/bioinfo/Nilesh/HRRdb_Samples/Scripts/columns36.csv") 
        
    #usecols to specify the columns to be read
    cancercol=collist['cancervar'][collist['cancervar'].notna()] #columns to be read in cancervar file
    cancervar= pd.read_csv(cancerloc, usecols=cancercol, sep='\t')
    
    vcfcol=collist['vcf_wo_art'][collist['vcf_wo_art'].notna()]
    vcfcol=[int(i) for i in vcfcol] #converting it into integers
    vcf= pd.read_csv(vcfloc, usecols=vcfcol, sep='\t')
    
    annocol=collist['multianno'][collist['multianno'].notna()]
    annovar= pd.read_csv(annoloc,usecols=annocol, sep='\t')

    # modifying vcf position values
    
    for i in range(len(vcf)):
        if len(vcf['REF'][i])>len(vcf['ALT'][i]):
            vcf['POS'][i]=vcf['POS'][i]+1
    
    print("VCF position values modified...")
    #rename annovar columns
    annovar=annovar.rename(columns={'Chr':'CHROM', 'Start': 'POS', 'Ref':'REF', 'Alt':'ALT','Gene.knownGene':'Ref.Gene' })
    
    #rename annovar columns
    cancervar=cancervar.rename(columns={'#Chr':'CHROM', 'Start': 'POS', 'Ref':'REF', 'Alt':'ALT','Gene.knownGene':'Ref.Gene' })
    
        
    cancervar['CHROM']=list(map(str, cancervar['CHROM']))
    cancervar['CHROM']='chr' + cancervar['CHROM']
    
    aicdf= pd.merge(annovar,cancervar, on= ['CHROM','POS','End','REF','ALT'])
    merged_df= pd.merge(vcf, aicdf, how='outer', on= ['CHROM','POS'])   
    merged_df.fillna('.', inplace = True)
     
    #preparing intervar_inhouse columns
    for i in range(31,59):
        merged_df[merged_df.columns[i]]=merged_df.columns[i]+ ":" + merged_df[merged_df.columns[i]]
    
    merged_df['intervar_inhouse']=merged_df[list(merged_df.columns[31:59])].apply(lambda x: ', '.join(x[x.notnull()]), axis = 1)

    ##Inserting columns 
    merged_df.insert(6, "IGV_link", value=None, allow_duplicates=False)
    merged_df.insert(7, "Mutant_allelic_burden_%", value=None, allow_duplicates=False)
    
    #add igv data in to the empty column
    merged_df['End'] =  merged_df['End'].replace('.', 0).fillna(0)
    merged_df['End'] = merged_df['End'].astype(int)
    merged_df['IGV_link'] = merged_df['CHROM'].astype(str) + [':'] + merged_df['POS'].astype(str) + ['-'] + merged_df['End'].astype(str)
    
    #add allele frequency data to empty column
    allele_freq=list(merged_df.columns[merged_df.columns.str.contains(':AF')])[0]
    
    if sample_type== 'DNA [Blood]': 
        merged_df[allele_freq]=merged_df[allele_freq].str.split(',', expand=True)[0].replace('.',0).fillna(0)
        merged_df['Mutant_allelic_burden_%'] = merged_df[allele_freq].astype(float)*100
        merged_df = merged_df.round({'Mutant_allelic_burden_%' : 0})
    else:
        merged_df[allele_freq] = merged_df[allele_freq].replace('.',0).fillna(0)
        merged_df['Mutant_allelic_burden_%'] = merged_df[allele_freq].astype(float)*100
        merged_df = merged_df.round({'Mutant_allelic_burden_%' : 0})
    

    #function to split columns containing multiple gene names sep by ";"
    def splitDataFrameList(df,target_column,separator):
    
        def splitListToRows(row,row_accumulator,target_column,separator):
            split_row = row[target_column].split(separator)
            for s in split_row:
                new_row = row.to_dict()
                new_row[target_column] = s
                row_accumulator.append(new_row)
        new_rows = []
        df.apply(splitListToRows,axis=1,args = (new_rows,target_column,separator))
        new_df = pd.DataFrame(new_rows)
        return new_df
    
    #spliting the 'Ref.Gene' column
    merged_df = splitDataFrameList(merged_df, 'Ref.Gene', ';')
    
    #performing a left merge with canonical data file
    merged_df = pd.merge(merged_df, canonical, how="left", on="Ref.Gene")
    
    #removing duplicates
    merged_df = merged_df.drop_duplicates()
 
    if sample_type=="DNA [Blood]": #"DNA [FE,FFPE]":
            #splitting the AD column into two columns 
            alter_depth=list(merged_df.columns[merged_df.columns.str.contains(':AD')])[0]
            merged_df[alter_depth]= merged_df[alter_depth].fillna('.')
            merged_df[alter_depth]= merged_df[alter_depth].replace('.','.,.')
            merged_df[['Ref_Depth', 'Mutant_Depth']] = merged_df[alter_depth].str.split(",", expand=True).drop([2],axis=1)
    else:
        #splitting the AD column into two columns 
        alter_depth=list(merged_df.columns[merged_df.columns.str.contains(':AD')])[0]
        merged_df[['Ref_Depth', 'Mutant_Depth']] = merged_df[alter_depth].str.split(",", expand=True)

    merged_df=merged_df.dropna(axis='columns', how='all')
    print(f + " : merged")
        
    #re-arranging the index
    cols=list(merged_df.columns)
    colind= list(collist['reindex_wo_art'][collist['reindex_wo_art'].notna()])
    colindex=list( [cols[int(i)] for i in colind] )
                        
    final_df=merged_df[colindex]           
    tot_var=len(final_df)
    
    #writing merged file
    output_path= dirpath + "/FE_merged/" + f + '_merged_output.csv'        
    final_df.to_csv(output_path, index=False)
    

####
time_taken = ("Time taken: {:.2f} seconds".format(time.time() - start_time))   
os.system("notify-send 'MultiAnno File Merger' 'Process Finished'")
msg = f"MultiAnno File Merger Finished @Nilesh_Mukherjee | {time_taken}"
print(msg) #for stdout and telegram-send --stdin
telegram_send.send(messages=[msg]) #for telegram notification
