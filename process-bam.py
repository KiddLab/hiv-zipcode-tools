import sys
import pysam

from optparse import OptionParser

###############################################################################
USAGE = """
python process-bam.py   --bam <path for bam file from insertion mapping>

Will output hit set to file in same dir.




"""
parser = OptionParser(USAGE)
parser.add_option('--bam',dest='bamFile', help = 'bam file')

(options, args) = parser.parse_args()



if options.bamFile is None:
    parser.error('bamFile not given')
###############################################################################


print options.bamFile

of = options.bamFile + '.align-parse.txt'
print 'writing output to',of
outFile = open(of,'w')


nl = ['#zipCode','readName','insChrom','insBP','insStrand','insMapQ','shearChrom','shearBP','insStrand','shearMapQ']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

bamFile = pysam.AlignmentFile(options.bamFile, "rb")
# get r1 and r2, primary...
r1 = None
r2 = None
for read in bamFile:
    if read.is_secondary is True:
        continue
    if read.is_read1 is True:
        if r1 == None:
            r1 = read
        else:
            print 'found r1 when have r1...???'
            print r1
            print read
            sys.exit()
    if read.is_read2 is True:
        if r2 == None:
            r2 = read
        else:
            print 'found r2 when have r2...???'
            print r2
            print read
            sys.exit()
    if r1 == None or r2 == None:
        continue # get next read
    else:
        # do the analysis
        name1 = r1.query_name
        name2 = r2.query_name
        if name1 != name2:
            print 'names do not match'
            print r1
            print r2
            sys.exit()
            
        # r1 should be to the genome
        r1Chrom = ''
        r1Pos = ''
        r1Dir = ''
        r2Chrom = ''
        r2Pos = '' # shear point
        r2Dir = ''
        
        qName = r1.query_name
        zc = qName.split(':')[0]
        
        if r1.is_unmapped is True:
            r1Chrom = 'NA'
            r1Pos = '.'
            r1Dir = '.'
            insPoint = '.'
            insOrient = '.'
            mapQ = -1
        else:
            r1Chrom = r1.reference_name
            refPos = r1.get_reference_positions()
            r1Start = refPos[0] +1 # will work in 1 based... 
            r1End = refPos[-1] + 1 # will work in 1 based... 
            mapQ = r1.mapping_quality
            if r1.is_reverse is True:
                r1Dir = '-'  # reverse strand, so in positive dir   
                insPoint = r1End
                insOrient = 'fwd'             
            else:
                r1Dir = '+'  # positive strand, so in reverse dir
                insPoint = r1Start
                insOrient = 'rev'             
        
        # now do to R2...
        if r2.is_unmapped is True:
            r2Chrom = 'NA'
            r2Pos = '.'
            r2Dir = '.'
            shearPoint = 0
            mapQR2 = 0
        else:
            r2Chrom = r2.reference_name
            refPos = r2.get_reference_positions()
            r2Start = refPos[0] +1 # will work in 1 based... 
            r2End = refPos[-1] + 1 # will work in 1 based... 
            mapQR2 = r1.mapping_quality
            
            if r2.is_reverse is True:
                r2Dir = '-'  # reverse strand, so in positive dir
                shearOrient = 'rev'
                shearPoint = r2End   
            else:
                r2Dir = '+'  # positive strand, so in reverse dir
                shearOrient = 'fwd'
                shearPoint = r2Start                 

            
#        print 'info!!!'        
#        print r1
#        print r1Chrom,r1Start,r1End,r1Dir,insPoint,insOrient,mapQ
        
#        print 'r2!!!!'
#        print r2
#        print r2Chrom,r2Start,r2End,r2Dir,shearPoint,r2Dir,mapQR2
        
        nl = [zc,qName,r1Chrom,insPoint,insOrient,mapQ,r2Chrom,shearPoint,shearOrient,mapQR2]
        nl = [str(j) for j in nl ]
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)
 
       ###### end do the analysis
        r1 = None
        r2 = None
 
        