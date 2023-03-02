#!/bin/bash
SECONDS=0

echo "PATHFINDER"

input="/home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/attempt2_192_Somatic.txt"
while IFS= read -r line
do
echo "Starting ${line}"

#%% Option1
#echo "For Line ${line}" >> ~/Pathfound_filtered3.txt
find /media/bioinfo/Basecare_s4/TMB_files/Annotation/research/TMB/ -type f -name "${line}*".vcf >> ~/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/attempt2_192_Somatic_Paths.txt    
# -printf '%h\n'     -not -path "*fe_*" | sort -u

#%% Option2
# readlink -f ${line} >> Pathfound.txt

echo "Ending ${line}"

done < "$input"


echo "################## ALL FILES ARE DONE ################"

#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"

notify-send "Pathfinder" "FilePath Gathering Finished"
telegram-send "Pathfinder finished attempt2_192_Somatic_Paths in $SECONDS seconds"


