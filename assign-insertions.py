import sys


from optparse import OptionParser

###############################################################################
USAGE = """
python assign-insertions.py   --collapsed <txt file from collapse-insertion-map.py>

Will output updated set to file in same dir.



"""
parser = OptionParser(USAGE)
parser.add_option('--collapsed',dest='collapsedFile', help = 'collapsed location file')


(options, args) = parser.parse_args()



if options.collapsedFile is None:
    parser.error('collapsedFile not given')

###############################################################################
def get_max(data):
    # assumes sorted
    if len(data) == 1:
        return 0
    if len(data) == 0:
        print 'none left, problem??'
        sys.exit()
    if data[0][4] > data[1][4]:
        return 0
    if data[0][4] == data[1][4]:
        if data[0][5] >= data[1][5]:
            return 0
    print 'unclear!'
    print data[0]
    print data[1]
    sys.exit()
###############################################################################
def to_del(zc,c,bp,insStrand,row):
    if row[0] == zc:
        return True
    if row[1] == c and row[2] == bp and row[3] == insStrand:
        return True
    return False
###############################################################################
print options.collapsedFile

ofn = options.collapsedFile + '.assigned'
print 'writing outout to',ofn
outFile = open(ofn,'w')
nl = ['#zipFamily','insChrom','insBP','insStrand','numShearPoints','totReads']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

data = []
inFile = open(options.collapsedFile,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    line[2] = int(line[2])
    line[4] = int(line[4])
    line[5] = int(line[5])
    
    data.append(line)
    
inFile.close()

print 'Read in initial set of %i lines' % len(data)

# sort

data.sort(key=lambda kv: kv[5],reverse=True)
data.sort(key=lambda kv: kv[4],reverse=True)

numWrote = 0

while True:
    if len(data) == 0:
        break
    maxIndex = get_max(data)
    zc = data[maxIndex][0]
    c = data[maxIndex][1]
    bp = data[maxIndex][2]
    insStrand = data[maxIndex][3]
    
    nl = data[maxIndex]
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    numWrote += 1
    
    didDel = True
    while didDel is True:
        didDel = False
        for i in range(len(data)):
            if to_del(zc,c,bp,insStrand,data[i]) is True:
                didDel = True
                del data[i]
                break
    print 'Out of del loop!'
    print len(data)
    
    print 'Wrote out',numWrote
    
        



outFile.close()


