# -*- coding: utf-8 -*-
from __future__ import division
import os
import math
'''
class File:
    def __init__(self,file_path):
        self.path = file_path

    def readInfo(self,sep):
'''
os.chdir('/home/owht/下载/data_from_web/temp/TEST')


def Igt(a):
    a = map(float,a.split())     # for a gene, calculate the Igt of them in each tissue
    sumA = sum(a)
    b = []
    for element in a:
        b.append(-(math.log((element/sumA),2)))
    return b

def getP(a):
    a = map(float,a.split())
    pSum = sum(a)                #for a tissue sample, calculate the PX distribution of its genes
    b = []                       
    for element in a:
        b.append(element/pSum)
    return b
     

def JSD(a,b):
    construct = [] 
    for i in xrange(0,len(a)):
        temp = (a[1]+b[1])*0.5
        construct.append(temp)
    aaa = KLD(a,construct)
    bbb = KLD(b,construct)
    return 0.5*aaa+0.5*bbb
    
def KLD(a,b):
    temp = 0
    if len(a)!=len(b):
        print "FALSE!! Error: the two inputs differ in length."
    else:
        for i in xrange(0,len(a)):
            p = math.log((a[i]/b[i]),2)
            temp += a[i]*p

    return temp
#--------------------------------------------------------------------------------------------------------------
#Calculate the I(g,t) for gene g in tissue t
# Input: (no rownames and no colnames)
#           tissue1   tissue2   tissue3 ......   tissue N 
#     gene1 E(g1,t1)  E(g1,t2)  E(g1,t3)         E(g1,tN)
#     gene2      ...............................   
#     gene3      .............................
#     gene4      .........................
#      ...       .........................
#     geneN      ............................... E(gN,tN)
#
# Output:(no rownames and no colnames)
#  .............. tissue N ...........
#     gene N ...  I(gN,tN) ... 

def getIgt(path):
    f = open(path,'r')
    fw = open("getIgt.result.txt","w")          
    line = f.readline()                         
    while line:                                 
        temp = map(str,Igt(line))               
        for i in xrange(0,len(temp)):                    
            fw.write(temp[i])                   
            fw.write(',')                       
        fw.write("\n")                          
        line = f.readline()                     
                                                
    f.close()                                   
    fw.close()
                                                
def getPFile(path):
    f = open(path,'r')
    fw = open("getPFile.result.txt","w")          
    line = f.readline()                         
    while line:                                 
        temp = map(str,getP(line))               
        for i in xrange(0,len(temp)):                    
            fw.write(temp[i])                   
            fw.write(',')                       
        fw.write("\n")                          
        line = f.readline()                     
                                                
    f.close()                                   
    fw.close()


'''
def main():
    profileA = map(int,(raw_input()).split())
    profileB = map(int,(raw_input()).split())

    print JSD(getP(profileA),getP(profileB))


if __name__ == '__main__':
    main()
'''
