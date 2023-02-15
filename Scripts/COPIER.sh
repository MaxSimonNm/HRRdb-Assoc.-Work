#!/bin/bash
SECONDS=0
echo "Copier"

input="/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/s4_Outside_VCFs_Paths.txt"
while IFS= read -r line
do
echo "Starting ${line}"

#mkdir ${line}

cp ${line} /home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/Total_VCFs

#rsync -avz -e "ssh -i ~/.ssh/id_rsa" basecare@10.187.28.53:${line} /home/bioinfo/Nilesh/HRRdb_Samples/HRR_G+_VCFs_Elucidata/

#scp -r basecare@10.187.28.53:${line} /home/bioinfo/Nilesh/HRRdb_Samples/HRR_G+_VCFs_Elucidata/

echo "Ending ${line}"
done < "$input"
echo "################## ALL FILES ARE DONE ################"
#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"
notify-send "Copier" "File Copying Finished"
telegram-send "Copier finished Somatic_VCFs in $SECONDS seconds"

