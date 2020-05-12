# usage:
# ./runSTAR.sh

# edit the range in {1..x} where x == number of transcriptome files
# e.g. 6 life stage files for Hermetia, so 6 is given as the argument
for i in {1..6} ;
do read1=$(echo USO-Hil-$i* | cut -d ' ' -f 1) &&
read2=$(echo USO-Hil-$i* | cut -d ' ' -f 2) &&

echo "***Working on sample USO-Hil-$i"***

# begin STAR call
star  --genomeDir ../STARindex/ \
--readFilesCommand gzcat \
--readFilesIn "$read1" "$read2" \
--outFileNamePrefix ../STARout/USO-Hil-"$i"_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
ulimit -n 10000 ; done
