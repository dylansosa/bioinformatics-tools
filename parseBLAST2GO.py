# @author dylan sosa
# 2/12/18

from Bio.Blast import NCBIXML
import sys, getopt
import numpy as np

def parseXML(blastOutput):
    result_handle = open(blastOutput)
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print 'Hit Accession Number: ', alignment.accession
                print 'Hit Definition: ', alignment.hit_def
                print 'Query Sequence: ',hsp.query
                print 'Subject Sequence: ',hsp.sbjct,'\n'

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
    parseXML(inputXML)
