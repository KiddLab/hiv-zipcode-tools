import sys


from optparse import OptionParser

###############################################################################
USAGE = """
python preprocess-insertion-shear-brkkpnts.py   --alignparsefam <txt file from add-family-to-align-parse.py>

Will output updated set to file in same dir.

Imposes '3 nucleotide rule' to consider unique shear points


"""
parser = OptionParser(USAGE)
parser.add_option('--alignparsefam',dest='alignParseFam', help = 'align parse file with family info')


(options, args) = parser.parse_args()



if options.alignParseFam is None:
    parser.error('alignParseFam not given')
    
###############################################################################
def check_to_skip(line):
# see if it is one that we want to skip
# no fam, no chrom, mismatch, etc
    if line[0] == 'NotAssigned':
        return True
    if line[3] == 'HIV-1_provirs':
        return True
    if line[7] == 'HIV-1_provirs':
        return True
    if line[7] == 'HIV-1_provirs':
        return True
    if line[5] != line[9]:
        return True  #diff dirs...
    if line[3] != line[7] != line[9]:
        return True #diff chroms
    if int(line[6]) <= 10:   # skip if unclear ins point
        return True 
###############################################################################

print options.alignParseFam

ofn = options.alignParseFam + '.3ntrule'
print 'puttign output to',ofn

ofstats = ofn + '.stats'
outStats = open(ofstats,'w')
dataRow = []
inFile = open(options.alignParseFam,'r')
outFile = open(ofn,'w')
for line in inFile:
    if line[0] == '#':
       header = line
       outFile.write(line)
       continue
    line = line.rstrip()
    line = line.split()
    toSkip = check_to_skip(line)
    if toSkip is True:
        continue
    if int(line[10]) <= 10:  # shear MapQ check
        continue


    line[8] = int(line[8]) # make shear bp int
    dataRow.append(line)
inFile.close()

print 'have %i data row to process' % len(dataRow)

# make list of uniqe chrom,pos,dir shear points
shearPosWorking = {}
for i in dataRow:
     c = i[7]
     p = i[8]
     d = i[9]
     if (c,p,d) not in shearPosWorking:
         shearPosWorking[(c,p,d)] = []

print 'setup %i unique c,p,d'  % len(shearPosWorking)

spdKeys = shearPosWorking.keys()
spdKeys.sort(key=lambda i: i[1])
spdKeys.sort(key=lambda i: i[0])


# ok, go through and get sets of everything that is within 3 of it -- then need to go through and a get a 

#setup for chaining -- do only self to self+3 (so just 1 side, taken advantage of fact that it is in sorted order

numMatch = 0
for i in range(len(spdKeys)):
    for j in range(i+1,len(spdKeys)):
        if spdKeys[i][0] != spdKeys[j][0]:
            continue
        if spdKeys[i][2] != spdKeys[j][2]:
            continue
            
        if  spdKeys[j][1] - spdKeys[i][1] <= 3:
            shearPosWorking[spdKeys[i]].append(spdKeys[j])
            numMatch += 1
        

print 'num match was ',numMatch

# setup new names
newName = {}
for i in range(len(spdKeys)):
    newName[spdKeys[i]] = '.'
print 'setup',len(newName)

# now, through and setup the chain...
tempCount = 0
for i in range(len(spdKeys)):
    if newName[spdKeys[i]] != '.':
        continue
    toSet = {}
    toSet[spdKeys[i]] = 1
    origLen = len(toSet)
    while True:
        k = toSet.keys()
        for j in k:
            match = shearPosWorking[j]
            for m in match:
                toSet[m] = 1
        newLen = len(toSet)
        if origLen == newLen:
            break
        origLen = len(toSet)
    # here -- we have the set finished
#    print 'toset is',toSet
    nl = []
    for j in toSet:
        newName[j] = spdKeys[i]
        t = j[0] + ':' + str(j[1]) + ':' + j[2]
        nl.append(t)
    nl = ','.join(nl)
    k = spdKeys[i][0] + ':' + str(spdKeys[i][1]) + ':' + spdKeys[i][2]
    ol = [k,str(len(toSet)),nl]
    ol = '\t'.join(ol) + '\n'
    outStats.write(ol)    
    tempCount += 1
    if tempCount % 500 == 0:
        print '\tdid',tempCount
    
for row in dataRow:
    c = row[7]
    p = row[8]
    d = row[9]
    
    n = newName[(c,p,d)]
    row[7] = n[0]
    row[8] = n[1]
    row[9] = n[2]
    row[8] = str(row[8])
    nl = '\t'.join(row) + '\n'
    outFile.write(nl)
            
outStats.close()
outFile.close()

print 'DONE'