import os, re
# wd = '/Users/dylansosa/Documents/SLU/Zhang/Winter-Break-paper_tree/ConvertToExcel/'
# os.chdir(wd)
# un-comment the above lines only if you are not using this script via BASH and want to test in an editor. 

# usage:
# python convert2csv.py > [output.csv]
with open('RHS.PID.toxins.domainarch', 'r') as data:
    for line in data:
        row = re.split(r'\s{2,}', line) # splitting if separated by >=2 whitespaces
        print row
# this output then needs to be opened in excel. Find and replace unnecessary characters.
