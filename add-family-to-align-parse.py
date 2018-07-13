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
zipFamDict = {}
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    zc = line[1]
    if zc == 'NotAssigned':
        continue
    zipFamList.append(zc)
    zipFamDict[zc] = 1
inFile.close()

print 'Read in %i zip families' % len(zipFamList)
print 'And the matching dict is %i long' % len(zipFamDict)  # to speed up lookups...

# add in two pass approach to make things go even faster.....

observedZipSet = {}
inFile = open(options.alignParse,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    zc = line.split()[0]
    observedZipSet[zc] = ''
inFile.close()
print 'Read in %i unique zips to match' % len(observedZipSet)

print 'Doing the preassignment...'
# now, go through each and assign zip to it
assigned = 0
notassigned = 0
n = 0
for zc in observedZipSet:
    did = False
    if zc in zipFamDict:
        did = True
        assigned += 1
        zipFam = zc
        observedZipSet[zc] = zipFam        
        continue
        
    
    for i in range(len(zipFamList)):
        # use lowmem=False to speed things up, since have reduced number of comparisons
        numMisMatches = zipcodetools.score_num_missmatches(zipFamList[i],zc)

        if numMisMatches <= options.editDistance:
            assigned += 1
            zipFam = zipFamList[i]
            did = True
            observedZipSet[zc] = zipFam            
            break            
    if did is False:
        notassigned += 1
        zipFam = 'NotAssigned'
        observedZipSet[zc] = zipFam            

print 'Assignment done:'
print 'Has assignment: %i  No assignment: %i' % (assigned,notassigned)        




numFound = 0
numNotFound = 0

ofn = options.alignParse + '.zipfamilies'
inFile = open(options.alignParse,'r')
print 'Matching to ',ofn
outFile = open(ofn,'w')
outFile.write('#%s\n' % options.zipFamilies)
for line in inFile:
    if line[0] == '#':
        line = line.rstrip()
        line = line[1:]
        line = '#zipFamily\t' + line + '\n'
        outFile.write(line)
        continue
    line = line.rstrip()
    zc = line.split()[0]
    zipFam = observedZipSet[zc]
    if zipFam == 'NotAssigned':
        numNotFound += 1
    else:
        numFound += 1
    
    nl = zipFam + '\t' + line + '\n'
    outFile.write(nl)

outFile.close()

print 'Num Matched in set %i' % numFound
print 'Num not matched in set %i' % numNotFound