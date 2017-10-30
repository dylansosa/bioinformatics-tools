### First version is for standard input
###

# @author dylan sosa
# Dr. Ahn
# BCB 5100
import os, sys, re

def reverse(s):
    """
    Reverse given DNA seq
    """
    return s[::-1]

def complement(s):
    """
    Return comeplementary DNA seq
    """
    complement = { "T": "A", "A": "T", "G": "C", "C": "G" }
    reversed = reverse(s)
    revcomplement = ''.join([complement[base] for base in reversed])
    print 'Reverse complement is:',revcomplement
    return revcomplement

def main():
    """
    Get input seq, check they are DNA nucleotides
    """
    seq = raw_input('Type DNA seq here: ').upper()
    nuc = ['A','T','G','C']
    for base in seq:
        if base not in nuc:
            print 'Please only enter DNA nucleotides'
            sys.exit()
    complement(seq)
    sys.exit()
    
if __name__ == '__main__':
    main()

###    
### Second version is for files
# @author dylan sosa
def reverseTranscript(s):
    complement = { "T": "A", "A": "T", "G": "C", "C": "G" }
    with open(s,'r') as inf:
        for line in inf:
            reverse = line[::-1]
            revcomplement = ''.join([complement[base] for base in reverse.strip().upper()])
        print revcomplement
        return revcomplement
    
reverseTranscript('file.txt')

###
### Third version is with Biopython
import os, sys, re
from Bio.Seq import Seq

def main():
    dna_seq = raw_input('Type DNA seq here: ')
    dna_seq = Seq(dna_seq)
    dna_seq = dna_seq.reverse_complement()
    print 'Reverse complement is: ', dna_seq
    sys.exit()

if __name__ == '__main__':
    main()

