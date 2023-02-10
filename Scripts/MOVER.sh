#!/bin/bash
SECONDS=0
echo "Mover"

input="Problem File List.txt"
while IFS= read -r line
do
echo "Starting ${line}"

#mkdir ${line}

mv /home/bioinfo/Nilesh/HRRdb_Samples/Multi_Anno_AllSamples/${line} /home/bioinfo/Nilesh/HRRdb_Samples/Problem_Samples/

echo "Ending ${line}"
done < "$input"
echo "################## ALL FILES ARE DONE ################"
#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"
notify-send "Mover" "File Moving Finished"
