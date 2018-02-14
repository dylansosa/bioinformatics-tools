# @author dylan sosa
# 2/12/18

from Bio.Blast import NCBIXML
import sys, getopt
import numpy as np

"""
What is happening here:
first returns the QUERY ID and SEQUENCE
then returns top 10 hit's accession numbers sorted by E-Value
"""

def parseXML(blastOutput):
    result_handle = open(blastOutput)
    blast_records = NCBIXML.parse(result_handle)
    count = 0
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print alignment.accession
                count += 1
                if count == 10:
                    exit(0)
                # print 'Hit Definition: ', alignment.hit_def
                # print 'E-Value', hsp.expect
                # print 'Query Sequence: ',hsp.query
                # print 'Subject Sequence: ',hsp.sbjct,'\n'

def getQueryID(blastOutput):
    result_handle = open(blastOutput)
    blast_records = NCBIXML.parse(result_handle)
    count = 0
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print '>'+record.query
                print hsp.query
                count += 1
                if count == 1:
                    return

def main(argv):
    """
    Parse XML format output from a local blast.
    """
    blastXML = ""
    try:
        opts, args = getopt.getopt(argv,"hx:",["blastXML="])
    except getopt.GetoptError:
        print '''
        heck
        '''
        sys.exit(2)

if __name__ == "__main__":
    arguments = main(sys.argv[1:])
    inputXML = sys.argv[2]
    getQueryID(inputXML)
    parseXML(inputXML)
