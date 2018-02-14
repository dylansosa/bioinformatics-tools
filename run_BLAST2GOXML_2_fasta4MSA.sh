#!/bin/bash

for f in ./blast2go_XML/*.xml ;

do python ./parseBLAST2GO_V2.py -x "${f}" > ./MSA/temp_bothID.txt ; # make ID file

queryID="$(cat ./MSA/temp_bothID.txt | egrep '^[>]' | tr -d '>')" ; # hold the query name as a variable in variable to be shared

cat ./MSA/temp_bothID.txt | head -2 > ./MSA/"$queryID"_tophits.fasta ; # make fasta file now so that query will be at top

cat ./MSA/temp_bothID.txt | sed -n '3,$p' > ./MSA/"$queryID"_hitsID.txt ; # get only hit ID to be used with e-utils

FS=$'\n'; for next in $(cat ./MSA/"$queryID"_hitsID.txt); do esearch -db protein -query $next | efetch -db protein -format fasta; done >> ./MSA/"$queryID"_tophits.fasta ; done
