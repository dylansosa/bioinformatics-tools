# -*- coding: utf-8 -*-
# !/usr/bin/env python2.7
# @author dylan sosa
# 2 November 2017
# BCB 5200
# Dr. Ahn
import sys, getopt
from Bio.Seq import Seq
from Bio import SeqIO

def main(argv):
    """
    Determines reverse complement of from the following file formats: .txt, .fasta, or input taken from user.
    """
    seqFile = ""
    method = ""
    output = ""
    try:
        opts, args = getopt.getopt(argv,"hf:m:o:",["seqFile=","method=","output="])
    except getopt.GetoptError:
        print '''\n-f: input seq file \n-m: the method to be used; default is biopython if file is given'''
        sys.exit(2)

    if len(sys.argv) == 1:
        print '\nNo argument given; result will be printed but not saved.'
        seq = raw_input('Type DNA seq here: ').upper()
        nuc = ['A','T','G','C']
        for base in seq:
            if base not in nuc:
                print 'Please only enter DNA nucleotides'
                sys.exit()
        complement(seq)
        sys.exit()

    for opt, arg in opts:
        if opt=="-h":
            print '''[Command line execution]\n-f <input seq file>\n-m enter either 'biopython' or 'personal' as your method\n-o specify output file'''
            sys.exit()
        elif opt in ("-f","--seqFile"):
            seqFile = arg
        elif opt in ("-m","--method"):
            method = arg
        elif opt in ("-o","--output"):
            method = arg

def reverse(f):
    """
    Reverse given seq
    """
    return f[::-1]

def complement(f):
    """
    Return comeplementary DNA seq
    """
    complement = { "T": "A", "A": "T", "G": "C", "C": "G" }
    reversed = reverse(f)
    revcomplement = ''.join([complement[base] for base in reversed])
    print 'Reverse complement is:',revcomplement
    return revcomplement

def biopythonReverseComplement(f):
    with open(f,'r') as dna_seq:
        seqList = []
        fasta = 0
        for line in dna_seq:
            if line.startswith('>'):
                ID = str(line.strip())
                fasta = 1
                continue
            seqList.append(line)
        seqString = ''.join(seqList)
        f_seqObject = Seq(seqString)
        f_reverse_seqObject = f_seqObject.reverse_complement()
        if fasta == 0:
            return str(f_reverse_seqObject[::-1]).strip()
        elif fasta == 1:
            revcomp = str(f_reverse_seqObject).strip()
            return ID + '(reverse)\n' + revcomp[::-1]

def reverseComplement(f):
    complement = { "T": "A", "A": "T", "G": "C", "C": "G" }
    revcomp = ''
    fasta = 0
    with open(f,'r') as inf:
        for line in inf:
            if line.startswith('>'):
                ID = str(line.strip())
                fasta = 1
                continue
                line = line[::-1]
            revcomplement = ''.join([complement[base] for base in line.strip().upper()])
            revcomp += '\n'+revcomplement
        if fasta == 0:
            return revcomp.strip()
        elif fasta == 1:
            return ID + '(reverse)\n' + revcomp.strip()

def writeOut(output,outputType):
    with open(output, 'w') as outfile:
        outfile.write(outputType)

if __name__ == "__main__":
    arguments = main(sys.argv[1:])
    if len(sys.argv) == 3: # if only file argument is given
        seqFile = sys.argv[2]
        method = 'biopython'
        output = '{}_reverse.txt'.format(sys.argv[2]) # output file name
        outputType = biopythonReverseComplement(seqFile) # determines which method's output is being used
        writeOut(output,outputType)
        sys.exit()

    elif len(sys.argv) == 5: # if given a file and method type
        seqFile = sys.argv[2]
        method = sys.argv[4]
        output = '{}_reverse.txt'.format(sys.argv[2]) # output file name
        outputType = biopythonReverseComplement(seqFile) # determines which method's output is being used
        writeOut(output,outputType)
        sys.exit()

    elif len(sys.argv) == 7: # if all three arguments are given
        seqFile = sys.argv[2]
        method = sys.argv[4]
        output = sys.argv[6] # output file name
        if method == 'biopython':
            outputType = biopythonReverseComplement(seqFile) # determines which method's output is being used
        elif method == 'personal':
            outputType = reverseComplement(seqFile)
        else:
            print '\n***Choose either biopython or personal as -m method choice.***\n'
        writeOut(output,outputType)
        sys.exit()
