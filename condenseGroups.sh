mkdir ./tempGrouping 
cp *.kalign tempGrouping/
cd tempGrouping
rm *.c.kalign
cat * | cut -d ' ' -f 1 | sort -u > ../c.gi

IFS=$'\n'; for line in $(cat ../c.gi) ; do egrep -rnwe "${line}" | cut -d ' ' -f 1 | cut -d ':' -f 1,2,3 >> compare ; done

cat compare| egrep -v "cat|compare|convert|fasta2seqrows|Option" > compare_clusters
cat compare_clusters | sort -u -t ':' -k3,3 | sort -u -t ':' -k1,1 > clusters2keep
# keep lines with unique ID
# then keep unique number of clusters 

cat clusters2keep | cut -d ':' -f 1 > files2use
sed -i '1d' files2use
# delete first line

mkdir ../groups2keep 
IFS=$'\n'; for line in $(cat files2use) ; do cp "${line}" ../groups2keep/ ; done
# move new files to work with to individual subdirectory 
