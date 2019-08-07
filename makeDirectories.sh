#!/bin/bash

for f in ./*.blastp ; do 
cDir="$(echo "${f}" | cut -d '.' -f 2 | sed 's/\///' )" ; # hold ID that began search for the C-term directory 
mkdir "$cDir" ;
mv "${f}" "$cDir" ;
csplit -f ./"$cDir"/cTerminal ./"$cDir"/"${f}" /Query=/ {*} ; 
rm ./"$cDir"/cTerminal00 ; done 

