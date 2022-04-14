#!/usr/bin/env python

import sys

with open("%s" %sys.argv[1],'r') as coord:
  file1=coord.readlines()

label="ATOM"
residue="VAL"
atom="CA"
chain="A"
count=1
outFile=open("%s" %sys.argv[2],'w')
for line in file1:
  line=line.strip()
  line=line.split()
  x=str(line[0])
  y=str(line[1])
  z=str(line[2])
  outFile.write(label.ljust(4))
  countStr=str(count)
  outFile.write(countStr.rjust(7))
  outFile.write(atom.rjust(4))
  outFile.write(residue.rjust(5))
  outFile.write(chain.rjust(2))
  outFile.write(countStr.rjust(4))
  outFile.write(x.rjust(12))
  outFile.write(y.rjust(8))
  outFile.write(z.rjust(8))
  outFile.write("\n");
  count+=1

outFile.close()
