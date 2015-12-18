
#计算CpG[O/E]#

from __future__ import division

thesequence = raw_input("Please input the sequence: ")

def numberOfCpG(sequence):
    return sequence.count("cg")

def numberOfC(sequence):
    return sequence.count("c")

def numberOfG(sequence):
    return sequence.count("g")

def CpGOE(sequence):
    NumberCpG = numberOfCpG(sequence)
    NumberOfC = numberOfC(sequence)
    NumberOfG = numberOfG(sequence)
    sequenceLength = len(sequence)
    GCContent = (NumberOfG+NumberOfC)/sequenceLength
    theCpgOE =  NumberCpG/(sequenceLength*(GCContent)**2)

    return theCpgOE
