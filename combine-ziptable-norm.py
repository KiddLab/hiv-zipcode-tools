#!/usr/bin/env python2
import sys
import gzip

from optparse import OptionParser

###############################################################################
USAGE = """
python combine-ziptable-norm.py  --out <path for output, must end in .ziptable.gz>
                                 --tables <names of *.ziptable.gz to combine > 

Assumes table is gziped, and that column 3 is (proper) normalized fraction.

Combines to new file with column 3 as normalized fraction, column 2 set to -1.



"""
parser = OptionParser(USAGE)
parser.add_option('--out',dest='outFile', help = 'merged and normed output file')
parser.add_option('--tables',dest='tables', help = 'files to merge')

(options, args) = parser.parse_args()



if options.outFile is None:
    parser.error('outFile not given')
if options.tables is None:
    parser.error('tables not given')
###############################################################################

toDo = []
toDo.append(options.tables)
for f in args:
    toDo.append(f)


if options.outFile[-12:] != '.ziptable.gz':
    print 'name error -- expect file %s to end in \'ziptable.gz\' ' % options.outFile
    sys.exit()

for f in toDo:
	if f[-12:] != '.ziptable.gz':
		print 'name error -- expect file %s to end in \'ziptable.gz\' ' % f
		sys.exit()
		
numFiles = len(toDo)
print 'will normalize for %i files' % numFiles

zipList = {}
for f in toDo:
    print f
    inFile = gzip.open(f,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        zc = line[0]
        fraction = float(line[2])
        if zc not in zipList:
            zipList[zc] = 0.0
        zipList[zc] += fraction            
    inFile.close()
    print 'Len:',len(zipList)

print 'Normalizing'
for k in zipList:
    zipList[k] = zipList[k] / float(numFiles)

zipKeys = zipList.keys()
zipKeys.sort(key= lambda k: zipList[k], reverse=True)    		

print 'Writing output to',options.outFile
outFile = gzip.open(options.outFile,'w')
for k in zipKeys:
    outFile.write('%s\t-1\t%.8f\n' % (k,zipList[k])) 
outFile.close()
print 'DONE'
