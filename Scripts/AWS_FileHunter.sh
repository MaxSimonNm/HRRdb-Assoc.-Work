#!/bin/bash
SECONDS=0

echo "
______  _  _         _   _                _               
|  ___|(_)| |       | | | |              | |              
| |_    _ | |  ___  | |_| | _   _  _ __  | |_   ___  _ __ 
|  _|  | || | / _ \ |  _  || | | || '_ \ | __| / _ \| '__|
| |    | || ||  __/ | | | || |_| || | | || |_ |  __/| |   
\_|    |_||_| \___| \_| |_/ \__,_||_| |_| \__| \___||_|   
                                                          
                                                          
"
                                                                     
echo "Hunting Files with maximum file size in AWS"

input="/home/bioinfo/Nilesh/HRRdb_Samples/Scripts/file.txt"
while IFS= read -r line
do
echo "Starting ${line}"

aws s3 ls s3://4basecare-databackup/ --profile=4basecare --recursive --human-readable --summarize | grep "${line}" | sed '/gz./d' | sort -k 3 -r | head -n 1 >> hunted.txt

echo "Ending ${line}"

done < "$input"

echo "################## ALL FILES ARE DONE ################"

#### Total time taken at this section
echo "Elapsed Time: $SECONDS seconds"

notify-send "AWS_FileHunter" "FileHunter Collected Files with maximum filesize"
telegram-send "AWS_FileHunter finished hunting down files with maximum filesize in $SECONDS seconds @Nilesh_Mukherjee"
