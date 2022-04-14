#!/usr/bin/env python

import sys
import subprocess
import os.path

longer_names={'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP',
              'C':'CYS', 'E':'GLU', 'Q':'GLN', 'G':'GLY',
              'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS',
              'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
              'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}

datadir=sys.argv[1]

ssDat=open("%s/seq.dat" %datadir,'w')
turnFile=open("%s/turn.txt" %datadir,'w')

#Parse input argument
inputFile=open("%s/input.txt" %datadir,'r')
inputLines=inputFile.readlines()
inputFile.close()
ssString = ""
count=1
for line in inputLines:
    line=line.strip() 
    aatype=line[0]
    sstype=line[2]
    ssString=ssString+sstype
    ssDat.write(str(count).rjust(5))
    ssDat.write(longer_names[line[0]].rjust(6))
    if(sstype=='H'):
        ssDat.write("2".rjust(5))
    elif(sstype=='E'):
        ssDat.write("4".rjust(5))
    else:
        ssDat.write("1".rjust(5))

    ssDat.write("9".rjust(5))
    ssDat.write("\n")
    count+=1
ssDat.close()

splitSS=list(str(ssString))
length=count-1
numTurns=0
probTurn=[]
for i in range(length):
    if(i<length-4):
        if(splitSS[i]=='C' and splitSS[i+1]=='C' and splitSS[i+2]=='C' and splitSS[i+3]=='C'):
            probTurn.append('0.333333')
            numTurns=numTurns+1
        else:
            probTurn.append('-0.333333')
    else:
        probTurn.append('-0.333333')

totcor=float(numTurns)/float(length)
for j in range(0,length+2):
    if(j==0):
        turnFile.write("%s\n" %length)
    elif(j<=length):
        index=j-1
        turnFile.write(str(index).rjust(6))
        turnFile.write('1'.rjust(3))
        turnFile.write(probTurn[index].rjust(10))
        turnFile.write("\n")
    else:
        turnFile.write("%.6f" %totcor)
        turnFile.write(str(numTurns).rjust(3))
        turnFile.write(str(length).rjust(3))
turnFile.close()
