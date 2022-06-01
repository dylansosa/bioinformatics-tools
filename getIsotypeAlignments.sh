for i in *.isotype.fa ; do
kalign -i $i -o isotypeAlignments/$i.kalign.fa ; done 
