import sys


from optparse import OptionParser

###############################################################################
USAGE = """
python combine-assignment-table.py   --assignfile <txt file of assignments>
                                     --alignparsefam <txt file from add-family-to-align-parse.py>
                                     --zipfamilies <file of zip families to match>

 
Will output updated set to file in same dir.



"""
parser = OptionParser(USAGE)
parser.add_option('--assignfile',dest='assignFile', help = 'assign file')
parser.add_option('--alignparsefam',dest='alignParseFam', help = 'alignparsefam file')
parser.add_option('--zipfamilies',dest='zipFamilies', help = 'zipFamilies file')

(options, args) = parser.parse_args()



if options.assignFile is None:
    parser.error('assignFile not given')
if options.alignParseFam is None:
    parser.error('alignParseFam not given')
if options.zipFamilies is None:
    parser.error('zipFamilies not given')
###############################################################################

print 'Read in zip info from', options.zipFamilies

inFile = open(options.zipFamilies,'r')
intSite = {}
totReads = {}

orderList = []
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    orderList.append(line)
    zc = line[1]
    intSite[zc]= '.\t.\t.\t.\t.'
    totReads[zc] = 0
inFile.close()

print 'Read in %i' % len(orderList)


# now count total reads.

inFile = open(options.alignParseFam,'r')
for line in inFile:
    if line[0] == '#':
         continue
    line = line.rstrip()
    line = line.split()
    zc = line[0]
    totReads[zc] += 1
inFile.close()

print 'Added reads'


inFile = open(options.assignFile,'r')
for line in inFile:
    if line[0] == '#':
         continue
    line = line.rstrip()
    line = line.split()
    zc = line[0]
    nl = [line[1],line[2],line[3],line[4],line[5]]
    nl = '\t'.join(nl)
    intSite[zc] = nl
inFile.close()

ofn = options.assignFile + '.statstable'
print 'writing to',ofn
outFile = open(ofn,'w')

nl = ['#originalRank','zipCode','TZ5','totReadsAll','insChrom','insBP','insStrand','numShearPoints','totReads']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

for i in orderList:
    zc = i[1]
    nl = i
    nl.append(totReads[zc])
    nl.append(intSite[zc])
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
outFile.close()

