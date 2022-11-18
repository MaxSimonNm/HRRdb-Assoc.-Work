#!/bin/bash
SECONDS=0

echo "PATHFINDER"


input=4lSam_list.txt
while IFS= read -r line
do
echo "Starting ${line}"

#option1
echo "For Line ${line}" >> ~/Pathfound.txt
find ./ -name "${line}*" -not -path '*fe_*'-printf '%h\n' | sort -u >> ~/Pathfound_filtered.txt

#option2
# readlink -f ${line} >> Pathfound.txt

echo "Ending ${line}"

done < "$input"


echo "################## ALL FILES ARE DONE ################"

#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"

notify-send "Pathfinder" "FilePath Gathering Finished"


