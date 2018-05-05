import os
import sys
import zipcodetools

from optparse import OptionParser

###############################################################################
USAGE = """
python select-clusters.py  --clustertable <cluster table output>
                           --stepsize <window size to consider, default = 1>
                           --cutoff <required fraction unique in window, default = 0.90>

Clustertable is output from form-clusters-from-table.py

Find first point where >= cutoff % of zipcodes in rank [i,  i+stepsize] merge into previous
zipcodes.  That is, where cutoff % are not unique.


"""
parser = OptionParser(USAGE)
parser.add_option('--clustertable',dest='clusterTable', help = 'cluster table')
parser.add_option('--stepsize',dest='stepSize', type='int',default=1,help = 'step size')
parser.add_option('--cutoff',dest='cutOff', type='float',default=0.90,help = 'cutoff')

(options, args) = parser.parse_args()

if options.clusterTable is None:
    parser.error('clusterTable not given')

###############################################################################

print 'input file',options.clusterTable
print 'step size',options.stepSize
print 'cutoff',options.cutOff

outFileName = options.clusterTable + '.sel.%i-%.2f.txt' % (options.stepSize,options.cutOff)



myData={}
myData['clusterTable'] = options.clusterTable
myData['stepSize'] = options.stepSize
myData['cutOff'] = options.cutOff

zipcodetools.select_clusters(myData)

print '%i selected zips written to %s' % (len(myData['selectedClusterList']),outFileName)
outFile = open(outFileName,'w')
for i in range(0,len(myData['selectedClusterList'])):
    outFile.write('%i\t%s\n' % (i+1,myData['selectedClusterList'][i]))

outFile.close()
