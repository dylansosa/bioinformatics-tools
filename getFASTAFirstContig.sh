for f in *.fna; do 


bedtools getfasta -fi $f -bed $f.fai.firstContig.gff -fo $f.firstContig.fa ; done
