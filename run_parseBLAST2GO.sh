#!/bin/bash 
    for f in ./blast2go_XML/*.xml; do python ./parseBLAST2GO.py -x "${f}" >> parsedBLAST2GOXML.txt; done
