for f in groups2keep/*.kalign ; do cat "${f}" | seqrows2fasta > "${f}".fa ; done

for g in groups2keep/*.fa ; do hmmscan /home/databases/pfam/Pfam-A.hmm "${g}" | hmmer2table -e .01 > "${g}".pfamList ; done

for h in groups2keep/*.fa ; do hmmscan /home/databases/sluprofiledb/HMMDB "${h}" | hmmer2table -e .01 > "${h}".sluList ; done 

for i in groups2keep/*List ; do cat "${i}" >> allList ; done

tfilter -i 1 -i 0 -i 2..5 allList | domain2architecture > all.arch 
