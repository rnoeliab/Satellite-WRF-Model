#! /bin/bash
# read the order files 
filename='/home/username/file_order.txt'

while read line; do 
echo $line
wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov$line" --header "Authorization: Bearer 10E74C60-F16A-11E8-B781-B0B1502E57AA" -P /home/username/output/
done< $filename 

echo "Finished the download"
