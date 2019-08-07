find . -name "cTerminal*" -exec sh -c '	for f ; do cat "$f" | blast2bounded | boundtable2fasta | kalign | fasta2seqrows > "$f".kalign |& tee > ./log.txt; done ' sh {} +

find . -name "cTerminal*.kalign" -exec sh -c ' for f ; do cat "$f" | seqrows2fasta | fasta2gi | gi2taxnode -complete | padtable > "$f".tax; done ' sh {} +

find . -type d -exec cp addTax.sh {} \; # get addTax in all subdirectories
for d in ./*/ ; do (cd "$d" && endRange=$(find -name "cTerminal*.kalign" | wc -l) && sed -i "s/[0-9][0-9]}/$endRange}/" addTax.sh && ./addTax.sh); done # Edits addTax to have range for each subdirectory and runs addTax. 
