import os
import sys
import zipcodetools
import gzip

from optparse import OptionParser

###############################################################################
USAGE = """
python quantify-clusters.py  --families <table of families to consider>
                             --table <processed table of zip counds to quantify (.ziptable.gz)>
                             --editdistance <edit distance to merge zipcodes, default=1)
                             --out <output file name for quantification results>
                           

Quantify observations for set of zipcodes

--families is output from select-clusters.py
Assumes --table file  is gziped, and that column 3 is (proper) normalized fraction.

"""
parser = OptionParser(USAGE)
parser.add_option('--families',dest='familiesFile', help = 'zipcode families file')
parser.add_option('--table',dest='zipTable', help = 'zipcode table table (.ziptable.gz)')
parser.add_option('--out',dest='outFile', help = 'output file')
parser.add_option('--editdistance',dest='editDistance', type='int',default=1,help = 'edit distance to merge zipcodes in clustering')

(options, args) = parser.parse_args()

if options.familiesFile is None:
    parser.error('familiesFile not given')
if options.zipTable is None:
    parser.error('zipTable not given')
if options.outFile is None:
    parser.error('outFile not given')


###############################################################################

if options.zipTable[-12:] != '.ziptable.gz':
    print 'name error -- expect file to end in \'ziptable.gz\' '
    sys.exit()
    
    
print 'edit distance set to',options.editDistance

outFile = open(options.outFile,'w')
outFile.write('#Edit distance %i\n' % options.editDistance)
outFile.write('#Zipcode families %s\n' % options.familiesFile)
outFile.write('#Zipcode counts %s\n' % options.zipTable)


# read in the zipcode families
familySet = []
inFile = open(options.familiesFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    familySet.append([line[0],line[1],0.0])    
inFile.close()
print 'Setup %i families to match\n' % len(familySet)
noAssignment = 0.0
totPerDepth = 0.0

inFile = gzip.open(options.zipTable,'r')
lineNum = 0
for line in inFile:
    line = line.rstrip()
    line = line.split()
    zipcode = line[0]
    fracDepth = float(line[2])
    totPerDepth += fracDepth
    did = False
    for i in range(len(familySet)):
        numMisMatches = zipcodetools.score_num_missmatches(familySet[i][1],zipcode)
        if numMisMatches <= options.editDistance:
            familySet[i][2] += fracDepth
            did = True
            break
    if did is False:
       noAssignment += fracDepth 
    lineNum += 1
    if lineNum % 50000 == 0:
        print '... Line number',lineNum
    
    
            
inFile.close()

for i in range(len(familySet)):
    outFile.write('%s\t%s\t%.8f\n' % (familySet[i][0],familySet[i][1],familySet[i][2]))
outFile.write('NotAssigned\tNotAssigned\t%.8f\n' % noAssignment)
outFile.close()
print 'Total depth encountered: %.8f' % totPerDepth
print 'Not assigned: %.8f' % noAssignment



