#!/usr/bin/env python2
import os
import sys
import zipcodetools
from optparse import OptionParser


###############################################################################
USAGE = """
python extract-zips.py  --fq1 <fq.gz file for read 1>  --fq2 <fq.gz file for read 2>
                        --outbase <base dir for output>   --name <name for sample set>
                        --LTRPrimer <inner primer of LTR sequence, default=GCCAATCAGGGAAGTAGCCTTGTGTGTGG>
                        --linkerPrimer <inner primer of linker, default=AGGGCTCCGCTTAAGGGAC>
                       
                         --lefttarget <left of zipcode, default=AGATCCAcatcgATggc >
                         --righttarget <right of zipcode, default=cAGATCAAGGATATCTT >
                         --ltrproviral <LTR to proviral sequence to filter, default=GTCTTCGTTGGGAGTGAATTAGCCCTTCCAGTCCCCCCTTTTCTTTTAAAAAGTGGC >
                                              
                       
                        
                        
Will  extract zipcde and info into outbase/name directory

Some parameters are hardcoded in.  Assumes various programs are in users path.                        


"""
parser = OptionParser(USAGE)
parser.add_option('--fq1',dest='fq1', help = 'fq.gz for read 1')
parser.add_option('--fq2',dest='fq2', help = 'fq.gz for read 2')
parser.add_option('--outbase',dest='outBase', help = 'base dir for output')
parser.add_option('--name',dest='name', help = 'name for this sample')

parser.add_option('--LTRPrimer',dest='LTRPrimer',default='GCCAATCAGGGAAGTAGCCTTGTGTGTGG',help = 'LTR primer')
parser.add_option('--linkerPrimer',dest='linkerPrimer',default='AGGGCTCCGCTTAAGGGAC',help = 'region right of zip to align')

parser.add_option('--lefttarget',dest='leftTarget',default='AGATCCAcatcgATggc',help = 'LTR primer')
parser.add_option('--righttarget',dest='rightTarget',default='cAGATCAAGGATATCTT',help = 'LTR primer')
parser.add_option('--ltrproviral',dest='LTRproviral',default='GTCTTCGTTGGGAGTGAATTAGCCCTTCCAGTCCCCCCTTTTCTTTTAAAAAGTGGC',help = 'LTR to proviral sequence')



(options, args) = parser.parse_args()

if options.fq1 is None:
    parser.error('fq1 not given')
if options.fq2 is None:
    parser.error('fq2 not given')
if options.outBase is None:
    parser.error('output base dir not given')
if options.name is None:
    parser.error('name not given')


###############################################################################

print 'Starting extraction for name',options.name

myData = {}
myData['fq1'] = options.fq1
myData['fq2'] = options.fq2
myData['outBase'] = options.outBase
myData['name'] = options.name
myData['LTRPrimer'] = options.LTRPrimer
myData['linkerPrimer'] = options.linkerPrimer
myData['LTRproviral'] = options.LTRproviral


myData['LTRPrimer'] = myData['LTRPrimer'].upper()
myData['linkerPrimer'] = myData['linkerPrimer'].upper()
myData['LTRproviral'] = myData['LTRproviral'].upper()




myData['minLTRPrimerMatch'] = len(myData['LTRPrimer']) - 3
myData['minlinkerPrimerMatch'] = len(myData['linkerPrimer']) - 3


myData['minLTRproviralMatch'] = len(myData['LTRproviral']) - 4


myData['leftTarget'] = options.leftTarget
myData['rightTarget'] = options.rightTarget

myData['leftTarget'] = myData['leftTarget'].upper()
myData['rightTarget'] = myData['rightTarget'].upper()

myData['minLeftMatch'] = len(myData['leftTarget']) - 2  # 2 since reads not overlapping..
myData['minRightMatch'] = len(myData['rightTarget']) - 2 
myData['minZipLen'] =  20 - 3
myData['maxZipLen'] =  20 + 3

zipcodetools.set_default_prog_paths(myData)
zipcodetools.setup_output_dir(myData)


zipcodetools.order_integration_fq(myData)


zipcodetools.print_integration_extraction_stats(myData)
print 'Stats output to file:',myData['integrationExtractStatsFile']

