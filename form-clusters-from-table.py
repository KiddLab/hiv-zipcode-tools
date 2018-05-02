import os
import sys
import zipcodetools

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



from optparse import OptionParser

###############################################################################
USAGE = """
python form-clusters-from-table.py  --table <processed table to cluster>
                                    --name <name for sample set>
                                    --maxrank <maximum number to consider, default=1000)
                                    --editdistance <edit distance to merge zipcodes, default=1)


Assumes table is gziped, and that column 3 is (proper) normalized fraction.
Makes some assumptions about cutoffs

"""
parser = OptionParser(USAGE)
parser.add_option('--table',dest='tableFile', help = 'table.gz for clustering')
parser.add_option('--name',dest='name', help = 'name for this sample')
parser.add_option('--maxrank',dest='maxRank', type='int',default=1000,help = 'maximum number to consider clustering')
parser.add_option('--editdistance',dest='editDistance', type='int',default=1,help = 'edit distance to merge zipcodes in clustering')


(options, args) = parser.parse_args()

if options.tableFile is None:
    parser.error('tableFile not given')
if options.name is None:
    parser.error('name not given')


###############################################################################

print 'Starting clustering for ',options.tableFile
print 'max to do',options.maxRank
maxMismatches = options.editDistance
print 'max edit distance',maxMismatches


if options.tableFile[-12:] != '.ziptable.gz':
    print 'name error -- expect file to end in \'ziptable.gz\' '
    sys.exit()
    
baseOutFileName =    options.tableFile[0:-12]
print 'base out file name is', baseOutFileName


myData = {}
myData['zipTable'] = options.tableFile
zipcodetools.read_ziptable_to_list(myData)



zipRank = []
zipCulmFractionReads = []
zipNumClusters = []
zipClusterSeqs = []

zipRank.append(0)
zipCulmFractionReads.append(0.0)
zipNumClusters.append(0)


# start with first zipcode
i = 0
zipRank.append(i+1)
zipCulmFractionReads.append(myData['zipList'][i][1] + zipCulmFractionReads[i] )
zipClusterSeqs.append([myData['zipList'][i][0]])
zipNumClusters.append(len(zipClusterSeqs))


print 'ziprank',zipRank
print 'culmFrac',zipCulmFractionReads
print 'numclusters',zipNumClusters
print 'seqs',zipClusterSeqs


for i in range(1,options.maxRank):
    if i >= len(myData['zipList']):  # in case we go over...
        break
    
    indexHit = -1 # index of prev seqs that it matches to form cluster of
    for j in range(len(zipClusterSeqs)):
        s2 = zipClusterSeqs[j][0]
        numMisMatches = zipcodetools.score_num_missmatches(myData['zipList'][i][0],s2)
        if numMisMatches <= maxMismatches:
            indexHit = j
            break
    
    # here, checked them all
    if indexHit == -1 : # did not have any hits
        zipRank.append(i+1)
        zipCulmFractionReads.append(myData['zipList'][i][1] + zipCulmFractionReads[i] )
        zipClusterSeqs.append([myData['zipList'][i][0]])
        zipNumClusters.append(len(zipClusterSeqs))
    else:  # had a merge!
        print 'MERGE!',i
        zipRank.append(i+1)
        zipCulmFractionReads.append(myData['zipList'][i][1] + zipCulmFractionReads[i] )
        zipClusterSeqs[indexHit].append(myData['zipList'][i][0])
        zipNumClusters.append(len(zipClusterSeqs))
        
    if i % 250 == 0:
        print '\tchecking rank %i ...' % i
        
        
print 'iteration',i
print 'ziprank',zipRank[-10:]
print 'culmFrac',zipCulmFractionReads[-10:]
print 'numclusters',zipNumClusters[-10:]
#print 'seqs',zipClusterSeqs
print '\n'

outTableFile = baseOutFileName + '.%s.%i.%i.out.txt' % (options.name,maxMismatches,options.maxRank)
print 'writing outinfo to',outTableFile    


print len(zipRank)
print len(zipCulmFractionReads)
print len(zipNumClusters)
print len(zipClusterSeqs)

outFile = open(outTableFile,'w')
for i in range(1,len(zipRank)):
    nl = [ zipRank[i] ]
    nl.append(zipNumClusters[i])
    nl.append(zipCulmFractionReads[i])
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)

outFile.write('\n#Clusters\n')
for i in range(len(zipClusterSeqs)):
    nl = [str(i+1)]
    seqs = zipClusterSeqs[i]
    nl.append(str(len(seqs)))
    seqs = ':'.join(seqs)
    nl.append(seqs)
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)    
outFile.close()


### plotting cmds #####
outPlotFile = baseOutFileName + '.%s.%i.%i.out.png' % (options.name,maxMismatches,options.maxRank)
print 'writing plot to',outPlotFile    

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(zipRank,zipNumClusters,'o',color='blue',markeredgecolor='blue')
ax.set_xlabel('Zipcode Rank')
ax.set_ylabel('Number of Zipcode Families',color='blue')
ax.set_title(options.name)
ax.set_xlim([0,len(zipRank)+5])
ax.set_ylim([0,len(zipRank)+5])
ax.tick_params('y', colors='blue')


ax2 = ax.twinx()
ax2.plot(zipRank,zipCulmFractionReads,'o',color='red',markeredgecolor='red')
ax2.set_ylabel('Cumulative Read Fraction', color='red')
ax2.tick_params('y', colors='red')
ax2.set_ylim([0,1.0])
ax2.set_xlim([0,len(zipRank)+5])


plt.savefig(outPlotFile)









