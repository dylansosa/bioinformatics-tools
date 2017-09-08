# @author dylan sosa
# charles hauser
# BINF 4410
# 5 April 2017
from Bio.Blast import NCBIXML
import sys, getopt
import numpy as np
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
    for opt, arg in opts:
        if opt=="-h":
            print '''
            heck
            '''
            sys.exit()
        elif opt in ("-x","--blastXML"):
            blastXML = arg

def parseXML(blastOutput):
    result_handle = open(blastOutput)
    blast_records = NCBIXML.parse(result_handle)
    #E_VALUE_THRESH = 0.04
    out = open('parsedXML.csv','w')
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                #if hsp.expect < E_VALUE_THRESH:
                print '****Alignment****'
                print 'Chlam ID:', record.query.split()[0]
                chlamID = record.query.split()[0]
                print 'ProtID in my DB:', alignment.title.split()[1]
                myProtID = alignment.title.split()[1]
                #print 'Sequence:', alignment.title
                print 'Length:', alignment.length
                if hsp.expect == 0:
                    hsp.expect == float("1e-315")
                print 'E-Value:', hsp.expect
                evalue = -1/np.log(hsp.expect)
                print 'Query  ',hsp.query[0:75] + '...'
                print 'Match  ',hsp.match[0:75] + '...'
                print 'Subject',hsp.sbjct[0:75] + '...','\n'
                out.write(str(chlamID)+','+str(myProtID)+','+str(evalue)+'\n')

if __name__ == "__main__":
    arguments = main(sys.argv[1:])
    inputXML = sys.argv[2]
    parseXML(inputXML)
