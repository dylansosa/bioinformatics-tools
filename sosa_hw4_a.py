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
