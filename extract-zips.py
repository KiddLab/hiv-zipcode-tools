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
parser.add_option('--target',dest='target', help = 'target region to align, flanking zipcode')


(options, args) = parser.parse_args()

if options.fq1 is None:
    parser.error('fq1 not given')
if options.fq2 is None:
    parser.error('fq2 not given')
if options.outBase is None:
    parser.error('output base dir not given')
if options.name is None:
    parser.error('name not given')
if options.target is None:
    parser.error('target fasta not given')


###############################################################################

print 'Starting extraction for name',options.name

myData = {}
myData['fq1'] = options.fq1
myData['fq2'] = options.fq2
myData['outBase'] = options.outBase
myData['name'] = options.name
myData['targetFA'] = options.target
myData['leftTarget'] = 'ACGAAGACAAGATATCCTTGATCTG'
myData['rightTarget'] = 'GCCATCGATGTGGATCTACCACACA'


zipcodetools.set_default_prog_paths(myData)
zipcodetools.setup_output_dir(myData)
zipcodetools.run_flash(myData)
zipcodetools.make_filtered_fasta(myData)
#zipcodetools.run_exonerate(myData)
zipcodetools.get_zipcode_noindel(myData)

