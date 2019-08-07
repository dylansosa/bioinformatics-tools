for f in ./*.full.fa ; do sed -i 's/-//g' "${f}" && cat "${f}" | splishrps -db /home/databases//rpsdb/allprofiles | blast2table -pcut .01 | tfilter -i 0..1 -i 6..7 -i 4 > "${f}".rpsList && hmmscan /home/databases/pfam/Pfam-A.hmm "${f}" | hmmer2table -e .01 | tfilter -i 0..4 > "${f}".pfamList && hmmscan /home/databases/sluprofiledb/HMMDB "${f}" | hmmer2table -e .01 | tfilter -i 0..4 > "${f}".sluList && for i in ./"${f}"*List ; do cat "${i}" >> "${f}".all ; done  ; done

for j in ./*.all ; do tfilter -i 1 -i 0 -i 2..5 "${j}" | domain2architecture | padtable | sort -k2,2 > "${j}".arch ; done
