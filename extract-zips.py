import os
import sys
import zipcodetools
from optparse import OptionParser


###############################################################################
USAGE = """
python extract-zips.py  --fq1 <fq.gz file for read 1>  --fq2 <fq.gz file for read 2>
                        --outbase <base dir for output>   --name <name for sample set>
                        --target <fasta of zipcode target>
                        
                        
Will  extract zipcde and info into outbase/name directory

Some parameters are hardcoded in.  Assumes various programs are in users path.                        


"""
parser = OptionParser(USAGE)
parser.add_option('--fq1',dest='fq1', help = 'fq.gz for read 1')
parser.add_option('--fq2',dest='fq2', help = 'fq.gz for read 2')
parser.add_option('--outbase',dest='outBase', help = 'base dir for output')
parser.add_option('--name',dest='name', help = 'name for this sample')


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
myData['leftTarget'] = 'ACGAAGACAAGATATCCTTGATCTG'
myData['rightTarget'] = 'GCCATCGATGTGGATCTACCACACA'

myData['minLeftMatch'] = len(myData['leftTarget']) - 1
myData['minRightMatch'] = len(myData['rightTarget']) - 1
myData['minZipLen'] =  20 - 3
myData['maxZipLen'] =  20 + 3

zipcodetools.set_default_prog_paths(myData)
zipcodetools.setup_output_dir(myData)
zipcodetools.run_flash(myData)
zipcodetools.read_flash_stats(myData)

zipcodetools.get_zipcode_noindel(myData)
zipcodetools.count_extracted_zips(myData)
zipcodetools.print_extraction_stats(myData)

print 'Stats output to file:',myData['extractStatsFile']
print 'Set of zipcodes output to',myData['zipTable']

