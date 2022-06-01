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

ID=$(echo $name | cut -f 1 -d '.')

spp=$(echo tRNA_gene_count_information/$ID*| rev | cut -f 4-  -d '.' | rev | cut -f 6 -d '_')

sed -i "s/^>/>$spp./" $name.$aa.isotype.fa ; done ;done
