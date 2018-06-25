import sys


from optparse import OptionParser

###############################################################################
USAGE = """
python add-family-to-align-parse.py   --alignparsefam <txt file from add-family-to-align-parse.py>

Will output updated set to file in same dir.



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

ofn = options.alignParseFam + '.collapsed'
print 'puttign output to',ofn
insData = {}

inFile = open(options.alignParseFam,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    toSkip = check_to_skip(line)
    if toSkip is True:
        continue
    
    zc = line[0]
    insC = line[3]
    insBP = line[4]
    insStrand = line[5]
    shearC = line[7]
    shearBP = line[8]
        
    if zc not in insData:
         insData[zc] = {}
    
    if (insC,insBP,insStrand) not in insData[zc]:
        insData[zc][(insC,insBP,insStrand)] = {}
        insData[zc][(insC,insBP,insStrand)]['reads'] = 0
    insData[zc][(insC,insBP,insStrand)]['reads'] += 1
    if (shearC,shearBP) not in insData[zc][(insC,insBP,insStrand)]:
         insData[zc][(insC,insBP,insStrand)][(shearC,shearBP)] = 0
    insData[zc][(insC,insBP,insStrand)][(shearC,shearBP)] += 1
        
    

inFile.close()

zcList = insData.keys()
zcList.sort()

outFile = open(ofn,'w')
nl = ['#zipFamily','insChrom','insBP','insStrand','numShearPoints','totReads']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

for zc in zcList:
    insPointList = insData[zc].keys()
    for insPoint in  insPointList:
        nl = [zc,insPoint[0],insPoint[1],insPoint[2],len(insData[zc][insPoint]) -1,insData[zc][insPoint]['reads']]
        nl = [str(j) for j in nl]        
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)

        
outFile.close()

