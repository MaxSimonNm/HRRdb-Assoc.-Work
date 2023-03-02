#!/bin/bash
SECONDS=0
echo "Mover"

input="ls.txt"
while IFS= read -r line
do
echo "Starting ${line}"

#mkdir ${line}

mv /home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/attempt2_192_Somatic_VCFs/${line}* /home/bioinfo/Nilesh/HRRdb_Samples/s4_HDD_TMB_VCFs/terminated_vcf/

echo "Ending ${line}"
done < "$input"
echo "################## ALL FILES ARE DONE ################"
#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"
notify-send "Mover" "File Moving Finished"
telegram-send "Mover finished in $SECONDS @Nilesh_Mukherjee"
