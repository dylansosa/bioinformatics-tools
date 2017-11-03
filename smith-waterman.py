# -*- coding: utf-8 -*-
# !/usr/bin/env python2.7
# @author dylan sosa
# 23 February 2017
# BINF 4410, Bioinformatics
# Dr. Charles Hauser

import sys, getopt
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
import numpy as np

def main(argv):
    """
    A python implementation of the Smith-Waterman algorithm for local sequence alignment.
    """
    sequence_1 = ""
    sequence_2 = ""
    scoring_matrix = ""
    gap_penalty_value = ""
    output_file = ""
    try:
        opts, args = getopt.getopt(argv,"hs1:s2:m:g:o:",["sequence_1=","sequence_2=","scoring_matrix=","gap_penalty_value","output="])
    except getopt.GetoptError:
        print '''
        [Command line execution]
        To execute the script from the command line, follow these guidelines:
        >smith-waterman.py –s1 seq1.fa –s2 seq2.fa –m blosum62 –g gap_penalty_value –o sw_out.txt

        Where:
        -s1, -s2: 2 sequences to be aligned
        -m: the scoring matrix employed;
        -g: gap penalty employed
        -o: output file
        '''
        sys.exit(2)
    for opt, arg in opts:
        if opt=="-h":
            print '''
            [Command line execution]
            To execute the script from the command line, follow these guidelines:
            >smith-waterman.py –s1 seq1.fa –s2 seq2.fa –m blosum62 –g gap_penalty_value –o sw_out.txt

            Where:
            -s1, -s2: 2 sequences to be aligned
            -m: the scoring matrix employed;
            -g: gap penalty employed
            -o: output file
            '''
            sys.exit()
        elif opt in ("-s1","--sequence_1"):
            sequence_1 = arg
        elif opt in ("-s2","--sequence_2"):
            sequence_2 = arg
        elif opt in ("-m","--scoring_matrix"):
            scoring_matrix = arg
        elif opt in ("-g","--gap_penalty_value"):
            gap_penalty_value = arg

def makeBlosum(scoring_matrix):
    """
    Based on user choice, make blosum dictionary.
    """
    scoring_matrix = scoring_matrix.upper()
    blosum = {}
    if scoring_matrix == ('BLOSUM50' or 'blosum50'):
        for key in MatrixInfo.blosum50:
            blosum[key] = MatrixInfo.blosum50[key]
            blosum[key[::-1]] = MatrixInfo.blosum50[key] #getting second half of matrix
        return blosum

    elif scoring_matrix == ('BLOSUM62' or 'blosum62'):
        for key in MatrixInfo.blosum62:
            blosum[key] = MatrixInfo.blosum62[key]
            blosum[key[::-1]] = MatrixInfo.blosum62[key] #getting second half of matrix
        return blosum

    else:
        raise Exception('Heck my guy, this is not a valid BLOSUM matrix!!! Call either 62 or 50.')

def aMatrixNamedAButAlsoOneNamedT(blosum):
    """
    Take chosen scoring matrix as parameter. Read in each of the two fasta files given by user.  Matrix A holds blosum score for each i,j based on the two files given. Fill A[i][j] with max(diag, upper, and left scores). Have the matrix T hold directions whence each score came.
    """
    with open(sequence_1,'r') as fasta_seq1:
        for line in fasta_seq1.readlines():
            line = line.strip()
            if line.startswith(">"):
                continue
            line.strip()
        fasta_seq1 = line
    with open(sequence_2,'r') as fasta_seq2:
        for line in fasta_seq2.readlines():
            line = line.strip()
            if line.startswith(">"):
                continue
            line.strip()
        fasta_seq2 = line

    m = len(fasta_seq1) # rows
    n = len(fasta_seq2) # cols

    A = np.zeros((m+1, n+1),dtype=int)
    T = np.chararray((m+1, n+1))
    T[:] = 'N'
    for i in range(1, m+1):
        for j in range(1, n+1):
            fasta_seq1_for_gettingscore = fasta_seq1[i-1]
            fasta_seq2_for_gettingscore = fasta_seq2[j-1]
            A[i][j] = blosum[(fasta_seq1_for_gettingscore,fasta_seq2_for_gettingscore)]
            blosum_score = A[i][j]

            # convert matrix that was filled with blosum scores
            diag_score = A[i-1][j-1] + blosum_score
            up_score = A[i-1][j] - gap_penalty_value
            left_score = A[i][j-1] - gap_penalty_value
            max_score = max(diag_score,up_score,left_score)
            if max_score <= 0:
                max_score = 0
            elif max_score == diag_score:
                T[i][j] = 'D'
            elif max_score == left_score:
                T[i][j] = 'L'
            elif max_score == up_score:
                T[i][j] = 'U'
            else:
                max_score = max_score
            # place the max score symbol in cell
            A[i][j] = max_score

    return A,T,m,n,fasta_seq1,fasta_seq2

def getMax(A,T,m,n,seq1,seq2):
        """
        Move through the matrix A and find the local maximum score. Return that and the coordinates.
        Conduct alignment of the given seqs.
        """
        local_max = 0
        local_max_coord = 0
        for i in range(1, m+1):
            for j in range(1, n+1):
                if A[i][j] < local_max:
                    local_max = local_max
                    local_max_coord = local_max_coord
                elif A[i][j] > local_max:
                    local_max = A[i][j]
                    local_max_coord = (i,j)
                elif A[i][j] == local_max:
                    local_max = A[i][j]
                    local_max_coord = (i,j)
        currentOptimalLocation = local_max_coord
        currentScore = local_max
        alignment1 = ''
        alignment2 = ''
        alignmentScore = 0
        G = np.chararray((m+1, n+1))
        G[:] = 'N'
        while currentScore != 0:
            if T[currentOptimalLocation] == 'U':
                alignment1 += seq1[currentOptimalLocation[0]-1]
                alignment2 += '-'
                currentOptimalLocation = (currentOptimalLocation[0]-1,currentOptimalLocation[1])
                G[currentOptimalLocation] = 'U'
            elif T[currentOptimalLocation] == 'D':
                alignment1 += seq1[currentOptimalLocation[0]-1]
                alignment2 += seq2[currentOptimalLocation[1]-1]
                currentOptimalLocation = (currentOptimalLocation[0]-1,currentOptimalLocation[1]-1)
                G[currentOptimalLocation] = 'D'
            elif T[currentOptimalLocation] == 'L':
                alignment1 += '-'
                alignment2 += seq2[currentOptimalLocation[1]-1]
                currentOptimalLocation = (currentOptimalLocation[0]-1,currentOptimalLocation[1])
                G[currentOptimalLocation] = 'L'
            alignmentScore += currentScore
            currentScore = A[currentOptimalLocation]
        rev_seq1 = alignment1[::-1]
        rev_seq2 = alignment2[::-1]

        print A
        print 'Local max is:',local_max,'at',local_max_coord,'\n\n','Aligned sequences:'
        print rev_seq1
        print rev_seq2
        print 'Alignment score:', alignmentScore
        return A,T,m,n,G,local_max_coord, alignmentScore, rev_seq1, rev_seq2

def prettify(A,T,m,n,G,LMC):
    """
    Use enumerated values to replace the locations in the matrix G with arrows and others!
    """
    UP_ARROW = u'\u21E7'
    LEFT_ARROW = u'\u21D0'
    DIAGONAL_ARROW = u'\u21D6'
    SPACER = u'\u25CB'
    BEGIN = u'\u2B50'
    symbolMatrix = []
    G[LMC] = 'B'
    for i in range(1, m+1):
        myString = ''
        for j in range(1, n+1):
            if G[i][j] == 'D':
                myString += DIAGONAL_ARROW+'  '
            elif G[i][j] == 'U':
                myString += UP_ARROW+'  '
            elif G[i][j] == 'L':
                myString += LEFT_ARROW+'  '
            elif G[i][j] == 'N':
                myString += SPACER+'  '
            else:
                myString += BEGIN+'  '
        symbolMatrix.append(myString)
    print '\nOptimal path:'
    for l in symbolMatrix:
        print l

def writeOut(out,score,RS1,RS2):
    with open(out, 'w') as outfile:
        outfile.write(RS1)
        outfile.write('\n')
        outfile.write(RS2)
        outfile.write('\n')
        outfile.write('Alignment Score is: ')
        outfile.write(str(score))

if __name__ == "__main__":
    arguments = main(sys.argv[1:])
    sequence_1 = sys.argv[2]
    sequence_2 = sys.argv[4]
    scoring_matrix = sys.argv[6]
    gap_penalty_value = int(sys.argv[8])
    output_file = sys.argv[10]
    print 'First FASTA file is: {0}'.format(sequence_1)
    print 'Second FASTA file is: {0}'.format(sequence_2)
    print 'Scoring matrix is: {0}'.format(scoring_matrix)
    A,T,m,n,seq1,seq2 = aMatrixNamedAButAlsoOneNamedT(makeBlosum(scoring_matrix))
    A,T,m,n,G,LMC,AS,RS1,RS2 = getMax(A,T,m,n,seq1,seq2)
    prettify(A,T,m,n,G,LMC)
    writeOut(output_file,AS,RS1,RS2)
