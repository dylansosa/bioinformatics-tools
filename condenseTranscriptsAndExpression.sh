#!/bin/bash
# need data file with header and one without. Use file without header as input
# adjust for loop range according to number of columns in data
# must have input file sorted on gene ID, not symbol

# command: condenseTranscripts2GeneExpression.sh
header="$(cat GSE99574_HiSAT2_dmel.transcript_level_TPM.FB.txt | head -n1)"

# if transcript expression values correspond to the same gene, sum and collapse into one row
for ((i = 4; i <=131; i++)); do awk -v column=$i '{ seen[$2] += $column } END { for (j in seen) print j, seen[j] }' test_Dmel_expression_noHeader | sort -k1 > Dmel_expression$i ; done

# gets just expression value column and saves in new files
for ((k = 4; k <=131; k++)); do awk '{print $2}' Dmel_expression$k > expression$k ; done
for n in $(seq 4 9); do mv expression$n expression0$n ; done;

# add expression columns to corresponding genes
cat Dmel_expression4 | cut -d ' ' -f 1 > geneList
paste -d '      ' geneList expression* > tmp_Dmel_expression.txt

# add and clean up header from second input file
sed "1 s/^/$header\n/" tmp_Dmel_expression.txt > header_Dmel_expr.txt
sed 's/SYMBOL   //g' header_Dmel_expr.txt > header_tmp_Dmel_expr.txt
sed 's/FBtrID   //g' header_tmp_Dmel_expr.txt > intermediate_Dmel_expression.txt
