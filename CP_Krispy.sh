echo "CP KRISPY"

input="Pathfound_cleaned.txt"
while IFS= read -r line
do
echo "Starting ${line}"

#mkdir -p /home/bioinfo/git/HRRdb/samples/${line}
cp -r /home/bioinfo/Nilesh/HRRdb_Samples/${line} /home/bioinfo/git/HRRdb/samples/

echo "Ending ${line}"
done < "$input"
echo "################## ALL FILES ARE DONE ################"
notify-send "Process Finished"


