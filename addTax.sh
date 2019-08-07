for i in {01..64} ;
do awk 'FNR==NR{a[NR]=$0;next}{$1=a[FNR]}1' cTerminal"${i}".kalign.tax cTerminal"${i}".kalign > ./cTerminal"${i}".c.kalign ; done
