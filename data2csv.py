import os, re

# usage:
# python convert2csv.py > [output.csv]
with open('RHS.PID.toxins.domainarch', 'r') as data:
    for line in data:
        row = re.split(r'\s{2,}', line) # splitting if separated by >=2 whitespaces
        print row
# this output then needs to be opened in excel. Find and replace unnecessary characters.
