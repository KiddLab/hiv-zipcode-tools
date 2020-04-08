#!/usr/bin/env python2
import os
import sys
import gzip

from optparse import OptionParser

###############################################################################
USAGE = """
python merge-quant-tables.py  --out <output table for merged results>
                              --in < set of inputs to merge>


Assumes --in files are the output of quantify-clusters.py

"""
parser = OptionParser(USAGE)
parser.add_option('--out',dest='outFile', help = 'output file for merged results')
parser.add_option('--in',dest='inFiles', help = 'files to merge')


(options, args) = parser.parse_args()

if options.outFile is None:
    parser.error('outFile not given')
if options.inFiles is None:
    parser.error('inFiles not given')

###############################################################################

toDo = []
toDo.append(options.inFiles)
for f in args:
    toDo.append(f)
    
print 'have %i to merge' % len(toDo)    
    
nameList = []
zipListOrder = []
zipCounts = {}

zipsToDo = {}

# first, read in order from file 1
inFile = open(toDo[0],'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    zipListOrder.append([line[0],line[1]])
    zipsToDo[line[1]] = 1
inFile.close()

print 'Found %i zips' % len(zipListOrder)


for fn in toDo:
    print fn
    name = fn.split('/')[-1].split('.')[0]
    print 'name is',name
    nameList.append(name)
    inFile = open(fn,'r')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        z = line[1]
        freq = line[2]
        if z not in zipsToDo:
             print 'ERROR!  found a zip not in set'
             print name
             print fn
             print line
             sys.exit()
        zipCounts[(z,name)] = freq    
    inFile.close()

print 'Read in all!'    

outFile = open(options.outFile,'w')
nl = ['#originalRank','zipCode']
nl.extend(nameList)
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

for i in zipListOrder:
    nl = [i[0],i[1]]
    z = i[1]
    for n in nameList:
        nl.append(zipCounts[(z,n)])
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
outFile.close()

