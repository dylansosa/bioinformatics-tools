for k in ./cTerminal[0-9][0-9].kalign  ; do (cat "${k}" | awk -F' ' '!seen[$1]++' | seqrows2fasta | kalign | fasta2seqrows > tmp && mv tmp "${k}") ; done
for k in ./cTerminal[0-9][0-9].kalign.tax  ; do (cat "${k}" | awk -F' ' '!seen[$1]++' > tmp && mv tmp "${k}") ; done
./addTax.sh
