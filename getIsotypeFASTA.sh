for i in *.highConfidence.fa  ; do

name=$(ls $i | rev | cut -f 4-  -d '.' | rev)
# name=$(ls $i | rev | cut -f 2- -d '.' | rev)
#egrep -wvf $i.pseudogeneCoordinates $name.bed > $i.highConfidence.bed ; done

#stdin=$(</dev/stdin)
# user input, selected isotype
# change eventually to use all isotypes
# maybe a file to look up


#egrep $1 $i > $name.highConfidence.$1

cat isotypes.txt | while read -r aa ; do

cat $i | egrep $aa -A1 --no-group-separator > $name.$aa.isotype.fa

spp=$(echo $name | cut -f 3 -d '_')

sed -i "s/^>/>$spp./" $name.$aa.isotype.fa ; done ;done

#for f in *.highConfidence.[a-z]*[A-Z][a-z]{2} ; do

#for f in *.highConfidence.*[A-Z]* ; do

#kalign -i $f -o $f.kalign.fa ; done 

#kalign -i $name.highConfidence.$1 -o $name.highConfidence.$1.kalign.fa ;done
#egrep $1 $i ; done

# get for all isotypes in a loop
# so make files for esch isotype 
