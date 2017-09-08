#!/usr/bin/env python2.7
'''
Created on Nov 30, 2015
Author: J. Dylan Sosa
In collaboration with: Isavannah-chan
Dr. Charles Hauser
BINF 3325
'''

import sys, getopt
import numpy as np
#import pylab as pl
import gzip
import re
from Bio import ExPASy
from Bio import SwissProt
#import matplotlib.pyplot as plt
#import numpy as np
#import plotly.plotly as py
#import progressbar
#from time import sleep
#import fish
from tqdm import tqdm


def main(argv):
    '''
    """
###############################################################
UCP Project, Mus musculus
--------------------------------------------------------------
sosa_UCP.py:
    A project designed to read microarray and expression value data, and extract values with pvalues < 0.05. These will be saved into a dictionary and then used to lookup accession IDs that will later be used to search BioDBnet. Finally, creates a sorted dictionary with uniprot IDs, GO terms, and fold change.

How to call:
    python sosa_UCP.py [-h] -e <expressionfilename> -c <chipfile> -b <bioDBfile>

-e  <expressionfile>
-c  <chipfile>
-b  <bioDBfile>
--------------------------------------------------------------
##############################################################
"""
'''
    expressionfile = ''
    chipfile = ''
    bioDBfile = ''
    try:
        opts, args = getopt.getopt(argv,'he:c:b:',['expressionfile=','chipfile=', 'bioDBfile='])
    except getopt.GetoptError:
        print'''How to call:
            python sosa_UCP.py [-h] -e <expressionfilename> -c <chipfile> -b <bioDBfile>

        -e  <expressionfile>
        -c  <chipfile>
        -b  <bioDBfile>'''
        sys.exit(2)
    for opt, arg in opts:
        if opt=='-h':
            print'''
            How to call:
                python sosa_UCP.py [-h] -e <expressionfilename> -c <chipfile> -b <bioDBfile>

            -e  <expressionfile>
            -c  <chipfile>
            -b  <bioDBfile>
            '''
            sys.exit()
        elif opt in ('-e','--expressionfile'):
            expressionfile = arg
        elif opt in ('-c','--chipfile'):
            chipfile = arg
        elif opt in ('-b','--bioDBfile'):
            bioDBfile = arg
    print 'The expression file is %s' % expressionfile
    print 'The chip file is %s' % chipfile
    print 'The BioDB file is %s' % bioDBfile
    print "\nComputing..."
    return expressionfile, chipfile, bioDBfile

def read_expressionfile(expressionfile):
    readFile = gzip.open(expressionfile, 'r')
    data = readFile.readlines()[1:] #skip ID line
    expressiondict= {}
    for line in data:
        expressionvalues = line.split() #WT value pval value pval value pval K5-UCP3 value pval value pval value pval
        expressiondict = load_expressiondict(expressionvalues, expressiondict)
    readFile.close()
    return expressiondict

def load_expressiondict(expressionvalues, expressiondict): #ID wtval wtpval wtval wtpval mutval mutpval mutval mutpval mut val mutpval
    WT_val = avg_val(expressionvalues[2], expressionvalues[4], expressionvalues[6], expressionvalues[1], expressionvalues[3], expressionvalues[5])
    K5UCP3_val = avg_val(expressionvalues[8], expressionvalues[10],expressionvalues[12],expressionvalues[7],expressionvalues[9], expressionvalues[11])
    if WT_val != 0 and K5UCP3_val != 0:
        foldchange = calc_foldchange(WT_val, K5UCP3_val)
        if abs(foldchange)>=2:
            expressiondict[expressionvalues[0]]= foldchange
    return expressiondict

def avg_val(pval1, pval2, pval3, val1,val2,val3):
    avg_val=0
    val1 = float(val1)
    val2 = float(val2)
    val3 = float(val3)
    pval1 = float(pval1)
    pval2 = float(pval2)
    pval3 = float(pval3)
    #import pdb;pdb.set_trace()
    #begin checking pvalues. Must be < 0.05 to be used.
    if pval1 < 0.05 and pval2 < 0.05 and pval3 < 0.05:
        avg_val = (val1+val2+val3)/3
    elif pval1 < 0.05:
        if pval2 <0.05:
            avg_val = (val1+val2)/2
        elif pval3 < 0.05:
            avg_val = (val1 + val3)/2
        else:
            avg_val = val1
    elif pval2 < 0.05:
        if pval3 < 0.05:
            avg_val = (val2+val3)/2
        else:
            avg_val = val2
    elif pval3 < 0.05:
        avg_val = val3
    return avg_val


def calc_foldchange(WT_value,K5UCP3_value):
    # foldchange = log(2)/log(treated/control)
    num = np.log(2)
    den = np.log(K5UCP3_value/WT_value)
    foldchange = num/den
    return foldchange

def matchIDs(expressiondict,chipfile):
    readfile = open(chipfile,'r')
    data = readfile.readlines()
    #import pdb; pdb.set_trace()
    for line in data:
        #import pdb; pdb.set_trace()
        if line[0].isdigit(): #we have reached the point where we can extract data in the file
            values = line.split()
            id = values[0]
            if id in expressiondict:
                accession = values[5]
                if accession.startswith('NM'):
                    expressiondict[accession] = expressiondict[id] #NM1234: 7.645452
                    del expressiondict[id]
    #import pdb; pdb.set_trace()
    for key in expressiondict.keys(): #removes the keys that did not match to a uniprot ID
        if key[0].isdigit():
            del expressiondict[key]
    return expressiondict

def match_uniprot_IDs(expressiondict,bioDBfile):
    readfile = open(bioDBfile,'r')
    data = readfile.readlines()
    #import pdb; pdb.set_trace()
    for line in data:
        values = line.split()
        id = values[0]
        accession = values[len(values)-1]
        if id in expressiondict:
            expressiondict[accession] = expressiondict[id]
            del expressiondict[id]
    for key in expressiondict.keys():
        if key.startswith('NM'):
            del expressiondict[key]
    return expressiondict
        #print values

def cross_ref_accessions(accessions):
    GOdict = {}
    for id in accessions:
        try:
            handle = ExPASy.get_sprot_raw(id)
            record = SwissProt.read(handle)
            #print record.description
            #print record.cross_references
            #all cross-references associated with the record INCLUDING Gene Ontology!
            #parse cross_references
            for ex_db_data in record.cross_references:
                #print ex_db_data
                extdb,extid,extdesc = ex_db_data[:3]
                if extdb=="GO" and extdesc.startswith("C"):
                    #print extdb,extid,extdesc
                    #GOdict[id] = [extdb,extid,extdesc]
                    GOdict[id] = [extdb,extid,extdesc]
            #break
        except:
            #print "WARNING: Handle %s not found" %id
            pass
    return GOdict

def dict_for_table(GOdict, expressiondict):
    for id in expressiondict:
        if id in GOdict:
            GOdict[id].append(str(expressiondict[id]))
    return GOdict

#def sort_dict_for_table(dict_for_table):
#    dict_for_table = [(key, value) for value, key in dict_for_table([(value,key) for key, value in dict_for_table.items()])]
#    return dict_for_table

def please_write(dict):
    with open('sorted_final.txt', 'w') as final_table:
    #    final_table.writelines('{}:{}\n'.format(k,v) for k, v in tableDict.items())
        final_table.write('UniprotIds\tFold Change\tGO Terms\tLocations\n')
        for key, value in sorted(tableDict.items(), key=lambda e: e[1][2].lower()):
            final_table.write(key+"\t")
            final_table.write(value[3]+"\t")
            final_table.write(value[0]+"\t")
            final_table.write(value[1]+"\t")
            final_table.write(value[2]+"\n")
    final_table.close()
    #return tableDict
    print'''
    The report has been succesfuly saved as <sorted_final.txt>!
          .-'"""""'-.
        .'           `.
       /   .      .    \
      :                 :
      |                 |
      :   \        /    :
       \   `.____.'    /
        `.           .'
          `-._____.-'

'''

#def please_sort(tableDict):
    #sortedvalues = sortedsorted(tableDict.values())
    #for values in sortedvalues:
    #    if values in
#    sorted(tableDict.items(), key=lambda e: e[1][2])
#    for key, value in sorted(tableDict.items(), key=lambda e: e[1][2].lower()):
#         print key, value
    #for key, value in sorted(tableDict.items(), key=lambda (k, v): v[2]):
    #     print key, value
#def histo(foldchange):
#    foldchange_data = foldchange
#    plt.hist(foldchange_data)
#    plt.title("Wt and UCP Fold Change Histogram")
#    plt.xlabel("Value")
#    plt.ylabel("Frequency")

#fig = plt.gcf()

#plot_url = py.plot_mpl(fig, filename='ucp_histo')

if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    arguments = main(sys.argv[1:])
    expressionfile = arguments[0]
    chipfile = arguments[1]
    bioDBfile = arguments[2]
    expressiondict = read_expressionfile(expressionfile)
    #for key,value in sorted(expressiondict.items()):
    #        print(key,value)
    expressiondict = matchIDs(expressiondict,chipfile)
    #for key,value in sorted(expressiondict.items()):
    #        print(key,value)
    expressiondict = match_uniprot_IDs(expressiondict,bioDBfile)
    #for key,value in sorted(expressiondict.items()):
    #        print(key)
    cross_ref_accessions(expressiondict.keys())
    GOdict = cross_ref_accessions(expressiondict.keys())
    #print GOdict
    tableDict = dict_for_table(GOdict, expressiondict)
    #for key,value in sorted(tableDict.items()):
        #print(key,value)
    #final_table = open('final_table.txt','w')
    #final_table.write(tableDict)
    #final_table.close()
    ##please_write(tableDict)
    #print sorted(tableDict.values())
    #please_sort(tableDict)
    please_write(tableDict)
    #makehisto(tableDict.values()[3])
###Make Graph
#d = tableDict
#X = np.arange(len(d))
#pl.bar(X, d.values(), align='center', width=0.5)
#pl.xticks(X, d.keys())
#ymax = max(d.values()) + 1
#pl.ylim(0, ymax)
#pl.show()
