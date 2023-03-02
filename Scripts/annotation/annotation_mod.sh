#!/bin/bash
SECONDS=0

#Enter the directory path where VCF is
cd /home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/s4_Outside_1-17_VCFs/

#Enter the VCF list
input="/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/outside-1-17-vcf_list.txt"
while IFS= read -r line
do
mkdir ${line}
echo "Starting ${line}"
#copying.hard-filtered.vcf from basespace and unzipping
#echo "##########Starting copying.hard-filtered.vcf from basespace and unzipping#################"
#cp {{projectdir}}/AppResults/${line}/Files/${line}.hard-filtered.vcf.gz ./
#gunzip ${line}.hard-filtered.vcf.gz
#echo "##########Ending copying.hard-filtered.vcf from basespace and unzipping#################"

#vcf to table
echo "##########Starting.hard-filtered.vcf conversion to tab#################"
python3 /home/bioinfo/Programs/VCF-Simplify-master/VcfSimplify.py SimplifyVCF -toType table -inVCF "${line}.hard-filtered.vcf" -outFile "${line}/${line}.tab"

#Put dummy coloumn with constant value and create final tab file required in filter engine 
awk '$++NF=NR==1?"Overlaps_from_BAM":01' ${line}/${line}.tab | cut -f 1-25 | sed 's/ /\t/g' > ${line}/${line}_final.tab

echo "##########Ending.hard-filtered.vcf conversion to tab#################"
#annovar
echo "##########Starting annovar annotation#################"
perl /home/bioinfo/Programs/Annotation_db/convert2annovar.pl -format vcf4old "${line}.hard-filtered.vcf" > "${line}/${line}.avinput"
perl /home/bioinfo/Programs/Annotation_db/table_annovar.pl "${line}/${line}.avinput" /home/bioinfo/Programs/Annotation_db/humandb/ -buildver hg19 -out "${line}/${line}_out" -remove -protocol ensGene,avsnp150,clinvar_20190305,intervar_20170202,intervar_20180118,esp6500siv2_all,exac03,gnomad211_exome,1000g2015aug_all,1000g2015aug_SAS,cadd13gt20,dbnsfp35a,dbscsnv11,dbnsfp31a_interpro -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring .
echo "##########Ending annovar annotation#################"
#preparing config file
perl /home/bioinfo/Nilesh/HRRdb_Samples/Scripts/annotation/config.pl ${line}/${line}
echo "##########Starting cancervar annotation#################"
#running cancervar
python3 /home/bioinfo/Programs/Annotation_db/CancerVar.py -c config.ini
echo "##########Ending cancervar annotation#################"
rm config.ini
done < "$input"
echo "################## ALL FILES ARE DONE ###########################"

#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"
notify-send "Annotator" "File Annotation Finished"
telegram-send "VCF Annotation Finished in $SECONDS @Nilesh_Mukherjee"

