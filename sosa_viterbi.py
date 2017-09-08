# -*- coding: utf-8 -*-
# !/usr/bin/env python2.7
# @author dylan sosa
# Dr. Charles Hauser
import re
from Bio import SeqIO
import numpy as np
import sys, getopt

def main(argv):
    membraneFile = ""
    solubleFile = ""
    stateSeqFile = ""
    trainingSetFile = ""
    try:
        opts, args = getopt.getopt(argv,"hm:s:q:t:",["membraneFile=","solubleFile=","stateSeqFile=","trainingSetFile="])
    except getopt.GetoptError:
        print '''[Command line execution]
        To execute the script from the command line, follow these guidelines:
        >lastname_viterbi.py –s membrane.txt –m soluble.txt –q stateSeq.txt –t training_set.txt

        Example:
        python testhmm.py -m set160_membrane.txt -s set160_soluble.txt -q set160_stateSequences.txt -t set160_training_set.txt
        '''
        sys.exit(2)
    for opt, arg in opts:
        if opt=="-h":
            print '''
            [Command line execution]
            To execute the script from the command line, follow these guidelines:
            >lastname_viterbi.py –s membrane.txt –m soluble.txt –q stateSeq.txt –t training_set.txt
            '''
            sys.exit()
        elif opt in ("-m","--membraneFile"):
            membraneFile = arg
        elif opt in ("-s","--solubleFile"):
            solubleFile = arg
        elif opt in ("-q","--stateSeqFile"):
            stateSeqFile = arg
        elif opt in ("-t","--trainingSetFile"):
            trainingSetFile = arg

def viterbi():
    membraneAACount = {}
    solubleAACount= {}
    transitionCount = {'SS': 0, 'SM': 0, 'MS': 0, 'MM': 0}
    membraneAAprobability = {}
    solubleAAprobability = {}
    stateSeqs = {}
    scoreMatrix = []
    testSeq = 'MLYNMNYLVFSLYKVFRQCLGCWRFLGIVKSPNT'

    with open(membraneFile, "r") as membrane:
        membrane = membrane.read()
        membrane = membrane.replace('\n', '')
        lenMembrane = len(membrane)

    with open(solubleFile, "r")as soluble:
        soluble = soluble.read()
        soluble = soluble.replace('\n', '')
        lenSoluble = len(soluble)

    with open(stateSeqFile, "r") as stateSeq:
        state = None
        s = ''
        membraneList = []
        solubleList = []
        for line in stateSeq.readlines():
            line = line.replace('o', 'S')
            line = line.replace('i', 'S')
            if line[0] == '>':
                if state:
                    stateSeqs[state] = s
                state = line[1:]
                stateSeqs[state] = ''
                s = ''
            else:
                s = s + line
                S = s.count('S')
                M = s.count('M')
                membraneList.append(M)
                solubleList.append(S)
        if state:
            stateSeqs[state] = s

        total = sum(membraneList) + sum(solubleList)
        avg = total/sum(membraneList)

    with open(trainingSetFile, "r")as trainingSet:
        trainingSet = trainingSet.read()
        trainingSet = trainingSet.rstrip('\n')

    countAA(membrane, membraneAACount)
    countAA(soluble, solubleAACount)
    countTransitions(stateSeqs, transitionCount)
    print 'Transition Frequencies:', transitionCount
    calculateP(membraneAACount, membraneAAprobability, lenMembrane)
    calculateP(solubleAACount, solubleAAprobability, lenSoluble)
    calculateFrequency(solubleAAprobability, membraneAAprobability, transitionCount, testSeq, scoreMatrix)

def countAA(file, dictionary):
    for aa in file:
        if aa in dictionary:
            dictionary[aa] += 1
        else:
            dictionary[aa] = 1

def countTransitions(stateSeqDict, transitionDict):
    soluble_S_list = []
    membrane_S_list = []
    membrane_M_list = []
    soluble_M_list = []
    solubleList = []
    membraneList = []
    length = 0
    for key in stateSeqDict:
        stateSeq = stateSeqDict[key]
        next_state = [stateSeq[i:i + 2] for i in range(len(stateSeq) - 1)]
        SS_count = next_state.count('SS')
        soluble_S_list.append(SS_count)
        SM_count = next_state.count('SM')
        membrane_S_list.append(SM_count)
        MM_count = next_state.count('MM')
        membrane_M_list.append(MM_count)
        MS_count = next_state.count('MS')
        soluble_M_list.append(MS_count)

    SS = float(sum(soluble_S_list))
    SM = float(sum(membrane_S_list))
    MM = float(sum(membrane_M_list))
    MS = float(sum(soluble_M_list))
    length = SS + SM + MM + MS
    SS = np.log2(float(SS / length))
    SM = np.log2(float(SM / length))
    MM = np.log2(float(MM / length))
    MS = np.log2(float(MS / length))
    transitionDict['SS'] = SS
    transitionDict['SM'] = SM
    transitionDict['MM'] = MM
    transitionDict['MS'] = MS

def calculateP(countDict, probDict, length):
    for aa in countDict:
        count = float(countDict.get(aa))
        prob = float(count/length)
        prob_norm = np.log2(prob)
        if aa not in probDict:
            probDict[aa] = prob_norm
        else:
            probDict[aa] = (int(countDict.get(aa)) / length)

def getCurrentState(P_current_is_S, P_current_is_M):
    if P_current_is_S > P_current_is_M:
        return 'S'
    else:
        return 'M'

def calculateFrequency(soluble_frequency, membrane_frequency, state_change_frequency, testSeq, scoreMatrix):
    stateSeq = ''
    total_score = 0
    for i in range (0, len(testSeq)):
        aa = testSeq[i]
        prev = testSeq[i-1]
        P_soluble = soluble_frequency.get(aa)
        P_membrane = membrane_frequency.get(aa)

        P_previous_is_S = soluble_frequency.get(prev)
        P_previous_is_M = membrane_frequency.get(prev)

        P_current_is_S = P_soluble + max(P_previous_is_S + state_change_frequency.get('SS'), P_previous_is_M + state_change_frequency.get('MS'))
        P_current_is_M = P_membrane + max(P_previous_is_S + state_change_frequency.get('SM'), P_previous_is_M + state_change_frequency.get('MM'))

        stateSeq += getCurrentState(P_current_is_S, P_current_is_M)

        if P_current_is_S > P_current_is_M:
            total_score += P_current_is_S
        else:
            total_score += P_current_is_M

    print 'Test Sequence is:', testSeq
    print 'State Sequence is:', stateSeq
    print 'Score:', total_score

if __name__ == "__main__":
    arguments = main(sys.argv[1:])
    membraneFile = sys.argv[2]
    solubleFile = sys.argv[4]
    stateSeqFile = sys.argv[6]
    trainingSetFile = sys.argv[8]
    viterbi()
