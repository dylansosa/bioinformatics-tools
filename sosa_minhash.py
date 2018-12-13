###########################################
# @author Dylan Sosa
# Dr. Chambers
# Algorithms in Computational Biology FA18
# December 11, 2018
###########################################
import sys, getopt
import math as m
import numpy as np
import pandas as pd
import random

def jaccardSimilarity(x, y):
    """
    Calculate Jaccard similarity index
    J(X,Y) = |XintersectY| / |XunionY|
    """
    numerator = len(set(x).intersection(y))
    denominator = len(set(x).union(y))
    return float(numerator) / denominator

def hashFunction(a, b, x, c = 402653189):
    """
    Awesome prime for this function found at the following address:
    https://planetmath.org/goodhashtableprimes
    """
    return ((a*x)+b)%c

def expected():
    return signatureSize*jaccardSimilarity(A,B)

def generateHashFamily():
    permutations = []
    for i in range(signatureSize):
        hashFunctionValues = (random.randint(0, (2**32)-1),random.randint(0, (2**32)-1))
        permutations.append(hashFunctionValues)
    return permutations

def minhash():
    """
    Python implementation of the minhash algorithm to measure string similarity.
    Compares Jaccard and minhash similarities.
    """
    signatureVector1 = np.full((signatureSize,1), float("inf"))
    signatureVector2 = np.full((signatureSize,1), float("inf"))
    for i, hashValues in enumerate(generateHashFamily()):
        a,b = hashValues
        for x in A:
            currentHash1 = hashFunction(a,b,x)
            if currentHash1 < signatureVector1[i]:
                signatureVector1[i] = currentHash1
        for y in B:
            currentHash2 = hashFunction(a,b,y)
            if currentHash2 < signatureVector2[i]:
                signatureVector2[i] = currentHash2

    matrix = np.concatenate((signatureVector1,signatureVector2),axis=1)
    signatureMatrix = pd.DataFrame(matrix, columns = ['set 1', 'set 2'])
    matchingValues = signatureMatrix.shape[0]*2 - np.unique(signatureMatrix).shape[0]
    print(signatureMatrix)
    print("Jaccard similarity is :", jaccardSimilarity(A,B))
    print("E value of minhash similarity is (assuming no collisions):", expected())
    print("Number of signature components is:", signatureSize)
    print("Minhash similarity is :", matchingValues / signatureSize)

A = set(['pero','jo','vull','correr','mon','vull','coneixer','gent'])
B = set(['jo','vull','anar','amb','el','vent','vull','coneixer','altres','terres'])
A = [ hash(x) for x in A ]
B = [ hash(x) for x in B ]

# A = [ hash(x) for x in C ]
# B = [ hash(x) for x in D ]
# C = {31,3,22,6,15,11}
# D = {15,30,7,11,28,3,17}

signatureSize = 10
minhash()
