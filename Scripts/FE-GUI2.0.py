import os
import pandas as pd
import warnings
import re



GUIpath='/home/bioinfo4/GUI2.0_Main_28th_Dec/GUI2.0-main/'
dirpath= '/home/bioinfo4/Patient_Samples/7th_Feb_LC-110_Urgent_samples/Filter_engine_checking/annotation/HRR'
sample_type= 'DNA[FFPE]'
test= 'HRR'

dirpath= dirpath + '/annotation/'

#importing external files for filter engine
canonical = pd.read_excel(GUIpath+ "/filter/canonical.xlsx", sheet_name=0, mangle_dupe_cols=True, engine='openpyxl')

######### selecting gene list #############
genes= pd.read_csv(GUIpath+ "/filter/genelist.csv")
testgenes= list(genes[test].dropna())
testgenes=[g.upper() for g in testgenes]

#making FE_merged and FE_filtered folders in the destination dir
folders= os.listdir(dirpath)

if 'FE_merged' in folders:
    os.system('rm -r ' + dirpath + '/FE_merged')
    os.system('rm -r ' + dirpath + '/FE_filtered')
folders= os.listdir(dirpath)    
os.system("mkdir " + dirpath + "/FE_merged")
os.system("mkdir " + dirpath + "/FE_filtered")
warnings.filterwarnings("ignore")

#Removing default file names from the sample name list
default_files=['annotation','ca','panel', 'dragen39.sh','cutadaptlog.txt','FQlog.txt','MSI','CNV','run_cnv.sh','cutadaptlog', 'QC','FE', 'FE_filtered','FE_merged','tmbmerged','tmbfiltered','tmb','sample','config.pl','annotation'] #pre-existing file names to be removed from the sample name list

for s in default_files:
    if s in folders:
        folders.remove(s)

filtered_df= pd.DataFrame(columns=['samplename','total_var','after exonic', 'after synony','after t4', 'after benign', 'after cadd', 'after pop_freq','after gen'])

#processing every sample in folder one by one
for f in folders:
    num=folders.index(f)
    f_path= dirpath + "/" + f
    files= os.listdir(f_path)
    cancer= [obj for obj in files if 'cancervar.hg19_multianno.txt.cancervar' in obj]
    annovar= [obj for obj in files if '_out.hg19_multianno' in obj]
    vcf= [obj for obj in files if 'final.tab' in obj]
    print(f)

    #locations of different files
    cancerloc= dirpath + "/" + f + "/" + cancer[0]
    vcfloc= dirpath + "/" + f + "/" + vcf[0]
    annoloc= dirpath + "/" + f + "/" + annovar[0]

    
    #Detecting dragen 3.6 or 3.9 & reading columns file
    vcfdetect=pd.read_csv(vcfloc, sep="/t")
    if 'INFO:hotspot' in vcfdetect.columns:
        collist= pd.read_csv(GUIpath+ "/filter/columns39.csv")
    else:
        collist= pd.read_csv(GUIpath+ "/filter/columns36.csv")        

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
    

    if sample_type=="DNA [Blood]":
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
    
    ############## Filtration ##########
    
    #######$$$$$$ Filtering for Clinvar pathogenic variants 
    clinvar= final_df[final_df['clinvar: Clinvar '].str.contains('pathogen', case=False, regex=True)]
    
    ####### filtering Func.knownGene
    if len(final_df)>0:
        df= final_df[final_df['Func.ensGene'].str.contains('exonic|splicing', case=False, regex=True)]
        df= df[df['Func.ensGene'].str.contains('RNA')==False]
        afknowngene= len(df)   
    else:
        afknowngene= 0
        
    ##### filtering synonymous
    if len(df)>0:
        df= df[df['ExonicFunc.ensGene']!='synonymous SNV']
        df['ExonicFunc.ensGene']= [i.title().replace(" Snv", " SNV") for i in df['ExonicFunc.ensGene']]
        afsynony= len(df)  
    else:
        afsynony= 0
    #####filtering pop freq
    
    popfreqs=['esp6500siv2_all','ExAC_ALL','ExAC_SAS','AF','AF_sas','1000g2015aug_all','1000g2015aug_SAS']
    
    print(f + " filtering in progress..")
    
    if len(df)>0:
        for p in popfreqs:
            df[p]=df[p].replace('.',0).fillna(0) 
            df[p]=  df[p].astype(float).round(4) #change 26-oct 12:45
            df=df[df[p]< 0.01] 
    
    afpop=len(df)
    ########## Removing Tier_IV

    if len(df)>0:
    
        df[' CancerVar: CancerVar and Evidence ']= [x.upper() for x in df[' CancerVar: CancerVar and Evidence ']]
        df=df[~df[' CancerVar: CancerVar and Evidence '].str.contains('BENIGN')]
        
    aft4=len(df)

    ######### Removing benign from intervar
    if len(df)>0:
        df['InterVar_automated']= [x.upper() for x in df['InterVar_automated']]
        df=df[~df['InterVar_automated'].str.contains('BENIGN')]    
        
    afben=len(df)
    
    if len(df)>0:
        df['CADD13_PHRED']=df['CADD13_PHRED'].replace('.',15).fillna(15)
        df['CADD13_PHRED']=[float(i) for i in list(df['CADD13_PHRED'])]
        df=df[df['CADD13_PHRED']>=15]
    
    afcad=len(df)
    
    
    ###### Gene filtering
    if len(df)>0:
        df['Ref.Gene']= [x.upper() for x in df['Ref.Gene']]
        df2=df[df['Ref.Gene'].str.contains('|'.join(testgenes))]
        df3= pd.DataFrame()
        for g in range(len(df2)):
            if ";" in df2['Ref.Gene'].iloc[g]:
                gene= df2['Ref.Gene'].iloc[g].split(";")[0]
            else:
                gene= df2['Ref.Gene'].iloc[g]
                
                if gene in testgenes:
                    df3=df3.append(df2.iloc[g])
        if len(df3)>0:  
            df3=df3[colindex] 
            afgen=len(df3)
        else:
            afgen=0
    else:
        df3=pd.DataFrame()
        afgen=0

    ##$$$$$ COMBINING clinvar variants and filtered variants together $$$$$###
    if len(df3)>0:
        df3=df3.append(clinvar,sort=False)
    
    #modified 27-10-2021
        df3['ExonicFunc.ensGene'] = df3['ExonicFunc.ensGene'].str.upper()
        df3['InterVar_automated'] = df3['InterVar_automated'].str.upper()
        df3[' CancerVar: CancerVar and Evidence '] = df3[' CancerVar: CancerVar and Evidence '].str.upper()
        
    #redundant code for clinvar data
        for p in popfreqs:
            df3[p]=df3[p].replace('.',0).fillna(0) 
            df3[p]=  df3[p].astype(float).round(4)
        
    ##Generating the Genomic Alt column
    
        #Rows not containing Frameshift Deletion
        df_n= df3[~df3['ExonicFunc.ensGene'].str.contains('deletion|insertion')]
        df_n['POS_x'] =  df_n['POS_x'].replace('.', 0).fillna(0)
        df_n['POS_x'] = df_n['POS_x'].astype(int)
        df_n['Genomic Alteration'] = df_n['CHROM_x'].astype(str) + [':g.'] + df_n['POS_x'].astype(str) + df_n['REF_x'] + ['>'] + df_n['ALT_x']  
        
        #Rows containing Frameshift deletion
        df_d= df3[df3['ExonicFunc.ensGene'].str.contains('deletion|insertion')]
        df_d['End_x'] =  df_d['End_x'].replace('.', 0).fillna(0)
        df_d['End_x'] = df_d['End_x'].astype(int)
        df_d['POS_x'] =  df_d['POS_x'].replace('.', 0).fillna(0)
        df_d['POS_x'] = df_d['POS_x'].astype(int)
        df_d['Genomic Alteration'] = df_d['CHROM_x'].astype(str) + [':g.'] + df_d['POS_x'].astype(str) + "_" + df_d['End_x'].astype(str)
        
        #merging both dataframes
        df3 = pd.merge(df_n, df_d, how='outer')
        #changing coloumn
        col_list = list(df3.columns)
        col_list[-2], col_list[-1] = col_list[-1], col_list[-2]
        df3 = df3[col_list]
        #modifying the ['AAChange.ensGene'] column 
        
        print("Modifying AA change....")
        
        aa=df3['AAChange.ensGene']
    
        to_replace={'A':'Ala','R':'Arg','N':'Asn','D':'Asp','B':'Asx','C':'Cys','E':'Glu','Q':'Gln','Z':'Glx','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}
    
        for a in aa:
            if 'ENS' in a:
                #splitting the list with , only if 'ENS' pattern is present- to avoid UNKNOWNS
                alist=a.split(',') #making a list with all the transcripts
                
                for al in alist: #running loop for changing each transcript
                        pattern = ":c.(.*?):p."
                        codon = re.search(pattern, al)
                        if codon is not None:
                            codon = re.search(pattern, al).group(1)
                        else:
                            codon=" "
                        #changing the codon nomenclature
                        if ('_' or 'del' or 'dup' or 'ins' or 'inv') in codon:
                            n_al=al #n_al where nothing is changed
                        
                        else:
                            ncodon=codon[1:len(codon)-1]+codon[0]+'>'+codon[len(codon)-1]
                            n_al=al.replace(codon, ncodon) #n_al with the replaced values
                        
                        #changing the protein nomenclature
                        if len(al.split(':p.'))>1:
                            prot=al.split(':p.')[1]
                        else:
                            prot=" "
                        for key, value in to_replace.items():
                            nprot = prot.translate(str.maketrans(to_replace))
                        
                        n_al=n_al.replace(prot, nprot)
                        
                        #making a list of aachange
                        aa=list(map(lambda st: str.replace(st, al, n_al), aa))   
            
        #replacing the modified values in the original dataframe
        df3['AAChange.ensGene'] = aa  
        
        print("AAChange.ensGene column modified ")
    
        df3=df3.dropna(axis='columns', how='all')
        df3=df3.astype(str)
        df3.drop_duplicates(subset=None, keep="first", inplace=True)
        
        ###Converting characters Cases in ExonicFunc.ensGene column###
        df3["ExonicFunc.ensGene"] = df3["ExonicFunc.ensGene"].str.capitalize()
        df3["ExonicFunc.ensGene"].replace("snv", "SNV", regex=True, inplace=True)
        df3["ExonicFunc.ensGene"].replace("deletion", "Deletion", regex=True, inplace=True) 
        df3["ExonicFunc.ensGene"].replace("insertion", "Insertion", regex=True, inplace=True) 
        ###Replacing if the last character of each cell is 'X' to 'Ter' in column AAChange.ensGene###
        df3['AAChange.ensGene'].replace("X$", "Ter", regex=True, inplace=True)
        
        df3['End_x'].replace("\.0", "", regex=True, inplace=True)  #varaible is object
        df3['CADD13_PHRED'].replace("\.0", "", regex=True, inplace=True)   
    
        df3['End_x']=df3['End_x'].astype(str).astype(int)  #varaible becomes int64 for ease in calculation
        for i in range(len(df3)):
            if (df3['ExonicFunc.ensGene'].iloc[i]=="Frameshift Insertion" or df3['ExonicFunc.ensGene'].iloc[i]=="Nonframeshift Insertion"):
                df3['End_x'].iloc[i]=df3['End_x'].iloc[i]+1
        df3['End_x'] = df3['End_x'].astype(str)
    
        df3.drop_duplicates(subset=None, keep="first", inplace=True)
        output_path= dirpath + "/FE_filtered/" + f + '_FENG.xlsx'  
        df3.to_excel(output_path, index=False)
        
        print(f + " :filtered")
        ###################################
    
####################------------------------------------------------------------------------------------
##################  Adding clinical sig and role data to the FENG file   #################################
#####################-------------- by prabir saha----------------------------------------------------------------------
    if len(df3)>0:
        DB_path_vus = GUIpath+"/filter/Clinical_Significance_DB.xlsx"   # Pathogen/VUS database
        DB_role = GUIpath + "/filter/Role_of_Gene_DB.xlsx"   # TSG/Oncogtenic database
        
        #### adding pathogenic VUS #######
        df1 = df3.reset_index(drop=True)
        df2 = pd.read_excel(DB_path_vus, sheet_name=0, engine='openpyxl')
        df2.drop(df2.filter(regex="Unname"),axis=1, inplace=True)
    
        row_num =  len(df1.index)
        
        list1 = []   # Empty list to put output of if else condition
        
        # Mathematical Strategy - If Ref allele length greter than Alt allele length then it consider deletion
        #                        - If Ref allele length less than Alt allele length then it consider insertion
        #                        - If Ref allele and Alt allele equal length then it consider SNV
        
        for i in range(0,row_num):
            if len(df1["REF_x"][i]) == len(df1["ALT_x"][i]):
                test_SNV = df1['CHROM_x'].astype(str)[i]+'|'+df1['POS_x'].astype(str)[i]+'|'+df1['REF_x'][i]+'|'+df1['ALT_x'][i]
                list1.append(test_SNV)
                
            else:
                if len(df1["REF_x"][i]) < len(df1["ALT_x"][i]):
                    test_col = df1["REF_x"][i]
                    test_len = len(test_col)
                    test_col2 = df1["ALT_x"][i]
                    s = df1['CHROM_x'].astype(str)[i]+'|'+df1['POS_x'].astype(str)[i]+'_'+df1['End_x'].astype(str)[i]+'|'+'ins'+ test_col2[test_len:] + "|"
                    list1.append(s)
                    
                else:
                    if len(df1["REF_x"][i]) > len(df1["ALT_x"][i]):
                        test_col1 = df1["ALT_x"][i]
                        test_len1 = len(test_col1)
                        test_col2 = df1["REF_x"][i]
                        t = df1['CHROM_x'].astype(str)[i]+'|'+df1['POS_x'].astype(str)[i]+'_'+df1['End_x'].astype(str)[i]+'|'+'del'+ test_col2[test_len1:] + "|"
                        list1.append(t)
        
        
        df1["New_Format"]= list1  #list1 data put in dataframe
        df1["New_Format"] = df1['New_Format'].str.replace("chr","")   # remove chr in the dataframe from new column
        
        #df1.to_excel(outputdir + "VC2E_new_format.xlsx", index=False)
        df_m = pd.merge(df1, df2, on ='New_Format', how ='outer')  # data featch from database
        my_column = df_m.pop('Clin Sig')
        df_m.insert(19,my_column.name, my_column)     # fetching column put in dataframe 
        
        input_row_01= len(df1.index)
        
        df_new_01 = df_m
        df_new_01["Clin Sig"] =  df_m["Clin Sig"].fillna(".")  #put dot in empty cell
        
        outp_row_01= len(df_new_01.index)
    
        for i in range(input_row_01,outp_row_01):
            df_new_01.drop([i],axis=0, inplace=True)   # remove unneccessary data
        
        
        df_new_01.rename(columns = {'Clin Sig':'Clin_Sig_inhouse'}, inplace= True)
        df_new_01.drop(df_new_01.filter(regex="Unname"),axis=1, inplace=True)
        
        #################### Role - TSG/Oncogenic ########################
        
        df = df_new_01
        de = pd.read_excel(DB_role, sheet_name=0, mangle_dupe_cols=True, engine='openpyxl')
        
        input_row= len(df.index)
        #print (input_row)
        df_n = pd.merge(df, de, on ='Ref.Gene', how ='outer')  #data fetching from database
        my_column_01 = df_n.pop('Role')
        df_n.insert(18,my_column_01.name, my_column_01)
        
        df_new = df_n
        df_new["Role"] =  df_n["Role"].fillna(".")
        
        outp_row= len(df_new.index)
        
        for i in range(input_row,outp_row):
            df_new.drop([i],axis=0, inplace=True)
        
        df_new.drop('New_Format',axis=1, inplace=True)
        df_new.drop('Therap_list',axis=1, inplace=True)
        df_new.drop('Pathway',axis=1, inplace=True)
    
        #df_new is the output of the above chunk of code #
        print('Clinical sig and role added to FENG file for ' + f)
    
        #####################--------------------------------------------------------------#######################    
        ######################## Sorting variants based on various criteria ##########################
        ######################-----------------------------------------------------------------####################
        
        df_sort= df_new
        if "#" in df_sort[' CancerVar: CancerVar and Evidence '][0]:
            cancervarscore=df_sort[' CancerVar: CancerVar and Evidence '].str.split(':', expand=True)[1].str.split('#', expand=True)[0]
        else:
            cancervarscore=df_sort[' CancerVar: CancerVar and Evidence ']
        
        ####converting columns to float
        df_sort[' CancerVar: CancerVar and Evidence ']=cancervarscore.astype(float)
        
        #according to 4basecare patho/vus
        df1=df_sort.replace({'Clin_Sig_inhouse': {'UNCERTAIN SIGNIFICANCE':1, '.' :2, 'VUS':3, 'DRUG RESPONSE':5, 'RISK FACTOR':5,'DRUG RESPONSE/PATHOGENIC':20, 'LIKELY PATHOGENIC':10, 'PATHOGENIC; DRUG RESPONSE':20, 'PATHOGENIC':20 }})         
        
        #scoring clinvar and intervar inhouse
        df1=df1.replace({'InterVar_automated':{'.':2,'UNCERTAIN_SIGNIFICANCE':3,'LIKELY_BENIGN':1,'BENIGN':1, 'LIKELY_PATHOGENIC':4,'PATHOGENIC':5}})
        df1=df1.replace({'clinvar: Clinvar ':{'clinvar: UNK ':2, 'clinvar: not_provided ':2, 'clinvar: other ':2, 'clinvar: Uncertain_significance ':3,'clinvar: Likely_benign ':1,'clinvar: Likely_pathogenic ':6,'clinvar: Pathogenic/Likely_pathogenic ':6,'clinvar: Pathogenic ':6}})
        
        ####scoring clinvar############
        
        #conflicting
        clinvar= df1['clinvar: Clinvar '].astype(str)
        conflict=list(filter(lambda x:'Conflicting' in x, clinvar))
        for c in list(set(conflict)):  
            df1=df1.replace({'clinvar: Clinvar ': {c : 3}})
        #pathogenic
        clinvar= df1['clinvar: Clinvar '].astype(str)
        pathogenic=list(filter(lambda x:'athogenic' in x, clinvar))
        for p in list(set(pathogenic)):
            df1=df1.replace({'clinvar: Clinvar ': {p : 6}})
        #benign
        clinvar= df1['clinvar: Clinvar '].astype(str)
        benign=list(filter(lambda x:'enign' in x, clinvar))
        for b in list(set(benign)):
            df1=df1.replace({'clinvar: Clinvar ': {b : 1}})
        #drug response
        clinvar= df1['clinvar: Clinvar '].astype(str)
        drugres=list(filter(lambda x:'rug_response' in x, clinvar))
        for d in list(set(drugres)):
            df1=df1.replace({'clinvar: Clinvar ': {d : 4}})

        df1['sort_score']= df1['clinvar: Clinvar '] + df1[' CancerVar: CancerVar and Evidence ']+ df1['Clin_Sig_inhouse']+ df1['InterVar_automated']
        df_sort['sort_score']=df1['sort_score']
        df_sort=df_sort.sort_values(by=['sort_score'], ascending=False)
    
        output_path= dirpath + "/FE_filtered/" + f + '_FENG.xlsx' 
        df_sort.to_excel( str(output_path), index=False)
    
    ###################################
    ### making filtered csv
    to_append= [f,tot_var,afknowngene,afsynony,afpop,aft4,afben,afcad,afgen]
    dflen=len(filtered_df)
    filtered_df.loc[dflen]=to_append
    print(to_append)  
    print( "###" + str(num+1) + " out of " + str(len(folders)) + " files done")

filtered_df.to_csv(dirpath+"/"+"FE_filtered.csv")        
print("#############################")
print("############ DONE ###########")
print("#############################")


############################################################_END FILTER_ENGINE ############################################################
