find . -maxdepth 1 -type d -exec cp condenseGroups.sh {} \; # add helper functions to subdirectories 
find . -maxdepth 1 -type d -exec cp cleanAlignments.sh {} \; 
find . -maxdepth 1 -type d -exec cp getArch.sh {} \;
for d in ./*/ ; do (cd "$d" && ./cleanAlignments.sh && ./condenseGroups.sh && ./getArch.sh) ; done

