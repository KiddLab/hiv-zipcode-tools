import sys
import zipcodetools


from optparse import OptionParser

###############################################################################
USAGE = """
python add-family-to-align-parse.py   --alignparse <txt file from process-bam.py>
                                      --zipfamilies <file of zip families to match>
                                      --editdistance <edit distance to merge zipcodes, default=1)
                                      

Will output updated set to file in same dir.



"""
parser = OptionParser(USAGE)
parser.add_option('--alignparse',dest='alignParse', help = 'align parse file')
parser.add_option('--zipfamilies',dest='zipFamilies', help = 'zip families file')
parser.add_option('--editdistance',dest='editDistance', type='int',default=1,help = 'edit distance to merge zipcodes in clustering')


(options, args) = parser.parse_args()



if options.alignParse is None:
    parser.error('alignParse not given')
if options.zipFamilies is None:
    parser.error('zipFamilies not given')
    
    
###############################################################################

print 'reading in zip families from',options.zipFamilies

inFile = open(options.zipFamilies,'r')
zipFamList = []
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    zc = line[1]
    if zc == 'NotAssigned':
        continue
    zipFamList.append(zc)
inFile.close()

print 'Read in %i zip families' % len(zipFamList)

numFound = 0
numNotFound = 0

ofn = options.alignParse + '.zipfamilies'
inFile = open(options.alignParse,'r')
print 'Matching to ',ofn
outFile = open(ofn,'w')
outFile.write('#%s\n' % options.zipFamilies)
numLine = 0
for line in inFile:
    if line[0] == '#':
        line = line.rstrip()
        line = line[1:]
        line = '#zipFamily\t' + line + '\n'
        outFile.write(line)
        continue
    line = line.rstrip()
    zc = line.split()[0]
    
    did = False
    for i in range(len(zipFamList)):
        numMisMatches = zipcodetools.score_num_missmatches(zipFamList[i],zc)
        if numMisMatches <= options.editDistance:
            numFound += 1
            zipFam = zipFamList[i]
            did = True
            break
            
    if did is False:
        numNotFound += 1
        zipFam = 'NotAssigned'
    nl = zipFam + '\t' + line + '\n'
    outFile.write(nl)

    numLine += 1
    if numLine % 100 == 0:
        print '... did',numLine,zipFam
    
    
outFile.close()

print 'Num Matched in set %i' % numFound
print 'Num not matched in set %i' % numNotFound