#!/usr/bin/python
###################
#J. Dylan Sosa
#Dr. Hauser
#18 September 2015
#BINF 3325
###################

import sys, getopt, gzip

# Get the total number of args, and arg list passed to the script
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

def main(argv):
    print"""
###############################################################
Protein Molecular Weight and Amino Acid Composition Calculator
--------------------------------------------------------------
molecWeightAAComp.py:
    Reads the protein sequence from a FASTA file
    and prints the sequence and ID information,
    the molecular weight in kDa, and the percent
    composition for the twenty common amino acids.


How to call:
    python molecWeightAAComp.py [-h] -i <filename>

-h    print instructional message
-i    <inputfile>
--------------------------------------------------------------
##############################################################
"""
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
    print 'Input file is: ', inputfile
    return inputfile # send variable 'inputfile' to program


if __name__ == "__main__":

    fin = main(sys.argv[1:])  #store inputfile as variable, 'fin'
#option 1 for accessing contents of file
    protseq_str=""
    protseq = gzip.open(fin,'r')
    for line in protseq: #print each line
        print(line)
        protseq_str+=line
#list of mw for aa
aaweights = {
    'A' : 71.09,  # alanine
    'R' : 156.19, # arginine
    'D' : 114.11, # aspartic acid
    'N' : 115.09, # asparagine
    'C' : 103.15, # cysteine
    'E' : 129.12, # glutamic acid
    'Q' : 128.14, # glutamine
    'G' : 57.05,  # glycine
    'H' : 137.14, # histidine
    'I' : 113.16, # isoleucine
    'L' : 113.16, # leucine
    'K' : 128.17, # lysine
    'M' : 131.19, # methionine
    'F' : 147.18, # phenylalanine
    'P' : 97.12,  # proline
    'S' : 87.08,  # serine
    'T' : 101.11, # threonine
    'W' : 186.12, # tryptophan
    'Y' : 163.18, # tyrosine
    'V' : 99.14   # valine
    }
weight = 0
#loop to begin adding mw
for aa in protseq_str:
    if aa in aaweights:
        weight = weight + aaweights[aa]
print("The molecular weight of the prot seq in kDa is", round(weight/1000,2))


###Begin percent composition by counting number of times each aa appears in seq
protseqLength = len(protseq_str)
#print len(protseq_str)
ACount = protseq_str.count('A')
RCount = protseq_str.count('R')
DCount = protseq_str.count('D')
NCount = protseq_str.count('N')
CCount = protseq_str.count('C')
ECount = protseq_str.count('E')
QCount = protseq_str.count('Q')
GCount = protseq_str.count('G')
HCount = protseq_str.count('H')
ICount = protseq_str.count('I')
LCount = protseq_str.count('L')
KCount = protseq_str.count('K')
MCount = protseq_str.count('M')
FCount = protseq_str.count('F')
PCount = protseq_str.count('P')
SCount = protseq_str.count('S')
TCount = protseq_str.count('T')
WCount = protseq_str.count('W')
YCount = protseq_str.count('Y')
VCount = protseq_str.count('V')

#want to print percent composition for each
percent_A = round(float(ACount)/float(protseqLength)*100,2)
print("\n")
print("The percent composition of A in the seq is:", percent_A)
percent_R = round(float(RCount)/float(protseqLength)*100,2)
print("The percent composition of R in the seq is:", percent_R)
percent_D = round(float(DCount)/float(protseqLength)*100,2)
print("The percent composition of D in the seq is:", percent_D)
percent_N = round(float(NCount)/float(protseqLength)*100,2)
print("The percent composition of N in the seq is:", percent_N)
percent_C = round(float(CCount)/float(protseqLength)*100,2)
print("The percent composition of C in the seq is:", percent_C)
percent_E = round(float(ECount)/float(protseqLength)*100,2)
print("The percent composition of E in the seq is:", percent_E)
percent_Q = round(float(QCount)/float(protseqLength)*100,2)
print("The percent composition of Q in the seq is:", percent_Q)
percent_G = round(float(GCount)/float(protseqLength)*100,2)
print("The percent composition of G in the seq is:", percent_G)
percent_H = round(float(HCount)/float(protseqLength)*100,2)
print("The percent composition of H in the seq is:", percent_H)
percent_I = round(float(ICount)/float(protseqLength)*100,2)
print("The percent composition of I in the seq is:", percent_I)
percent_L = round(float(LCount)/float(protseqLength)*100,2)
print("The percent composition of L in the seq is:", percent_L)
percent_K = round(float(KCount)/float(protseqLength)*100,2)
print("The percent composition of K in the seq is:", percent_K)
percent_M = round(float(MCount)/float(protseqLength)*100,2)
print("The percent composition of M in the seq is:", percent_M)
percent_F = round(float(FCount)/float(protseqLength)*100,2)
print("The percent composition of F in the seq is:", percent_F)
percent_P = round(float(PCount)/float(protseqLength)*100,2)
print("The percent composition of P in the seq is:", percent_P)
percent_S = round(float(SCount)/float(protseqLength)*100,2)
print("The percent composition of S in the seq is:", percent_S)
percent_T = round(float(TCount)/float(protseqLength)*100,2)
print("The percent composition of T in the seq is:", percent_T)
percent_W = round(float(WCount)/float(protseqLength)*100,2)
print("The percent composition of W in the seq is:", percent_W)
percent_Y = round(float(YCount)/float(protseqLength)*100,2)
print("The percent composition of Y in the seq is:", percent_Y)
percent_V = round(float(VCount)/float(protseqLength)*100,2)
print("The percent composition of V in the seq is:", percent_V)


protseq.close() #close file after done
