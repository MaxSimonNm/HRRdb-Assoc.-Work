#!/bin/bash
SECONDS=0
echo "
______         _          _                 
| ___ \       | |        | |                
| |_/ / _   _ | |_   ___ | |__    ___  _ __ 
| ___ \| | | || __| / __|| '_ \  / _ \| '__|
| |_/ /| |_| || |_ | (__ | | | ||  __/| |   
\____/  \__,_| \__| \___||_| |_| \___||_|   "

input="/home/bioinfo/Nilesh/HRRdb_Samples/Scripts/hunted.txt"
while IFS= read -r line
do
echo "Starting ${line}"

echo ${line} | sed 's/^[^ ]* [^ ]* [^ ]* [^ ]* //'  >> ./AWS_Hunted_Duplicate_VCFs_ButcheredPath.txt

echo "Ending ${line}"
done < "$input"
echo "################## ALL FILES ARE DONE ################"
#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"
notify-send "AWS_Butcher" "Filepaths Butchered"
telegram-send "AWS_Butcher finished butchering Hunted Files's path in $SECONDS seconds @Nilesh_Mukherjee"
