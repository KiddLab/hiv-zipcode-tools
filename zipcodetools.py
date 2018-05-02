# zipcodetools.py
# set of functions for working with revised HIV zipcodes
# python 2.X code
import sys
import os
import gzip
import subprocess

#####################################################################
# check to see if program is in PATH
# copied from https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
# is easier in Python 3 apparently...
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
#####################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print 'command failed'
        print cmd
        sys.exit(1)
#####################################################################
##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this could be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
#####################################################################
# setup paths to default programs to use....
def set_default_prog_paths(myData):
    print 'Checking for required programs\n'

    
    print 'Checking flash...'
    if which('flash') is None:
        print 'flash not found in path! please fix (module load?)'
        sys.exit()

#    print 'Checking exonerate...'
#    if which('exonerate') is None:
#        print 'exonerate not found in path! please fix (module load?)'
#        sys.exit()
#####################################################################
def setup_output_dir(myData):
    if myData['outBase'][-1] != '/':
        myData['outBase'] += '/'
    myData['outDir'] =  myData['outBase'] + myData['name']
    if os.path.isdir(myData['outDir']) is True:
        print 'Out dir exists'
    else:
        cmd = 'mkdir ' + myData['outDir']
        print cmd
        runCMD(cmd)
    myData['outDir'] += '/'
#####################################################################
def run_flash(myData):
    myData['flashDir'] = myData['outDir'] + 'flash-out/'    
    myData['flash_1'] = myData['flashDir'] + 'out.notCombined_1.fastq.gz'
    myData['flash_2'] = myData['flashDir'] + 'out.notCombined_2.fastq.gz'
    myData['flash_Frag'] = myData['flashDir'] + 'out.extendedFrags.fastq.gz'
    myData['flash_stats'] =  myData['flashDir'] + 'flashStats.log'
    
    if os.path.isdir(myData['flashDir']) is False:  # dir needs to exist at start for tee to work
        print 'making flash output dir!'
        cmd = 'mkdir ' + myData['flashDir']
        runCMD(cmd)
    
    if os.path.isfile(myData['flash_1']) is False:
        cmd = 'flash'  +' -r 75 -f 70 -s 2 --cap-mismatch-quals -z -t 1 ' + ' %s %s -d %s | tee %s' % (myData['fq1'],myData['fq2'],myData['flashDir'],myData['flash_stats'])
        print cmd
        runCMD(cmd)
    else:
        print 'looks like flash already ran!'
#####################################################################
def read_flash_stats(myData):  # read in data from file
    inFile = open(myData['flash_stats'],'r')
    lines = inFile.readlines()
    inFile.close()
    lines = lines[-10:]
    lines = [i.rstrip() for i in lines]
    
    if 'combination' not in lines[0]:
        print 'error -- log file for flash not formatted as expected!'
        for i,j in enumerate(lines):
            print i,j
        sys.exit()
        
    myData['flashInfo'] = {}
    myData['flashInfo']['totalpairs'] = int(lines[1].split()[-1])
    myData['flashInfo']['combinedpairs'] = int(lines[2].split()[-1])
    myData['flashInfo']['uncombinedpairs'] = int(lines[3].split()[-1])
    myData['flashInfo']['fraccombined'] = float(myData['flashInfo']['combinedpairs']) / float(myData['flashInfo']['totalpairs'])            
#####################################################################


#####################################################################
# to read in fastq like record
def get_4l_record(myFile):
    #fastq style file...
    # just return sequence len
    # -1 if last record
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    return [myLine1,myLine2,myLine3,myLine4]
#####################################################################
def make_filtered_fasta(myData):
    myData['flash_Frag_fasta'] = myData['flashDir'] + 'out.extendedFrags.fasta'
    if os.path.isfile(myData['flash_Frag_fasta']) is True:
        print 'Looks like already ran filter fasta'
        return
        
    inFile = gzip.open(myData['flash_Frag'],'r')
    outFile = open(myData['flash_Frag_fasta'],'w')
    n = 0
    print 'Filtering fastq...'
    while True:
        rec = get_4l_record(inFile)
        if rec == '':
            break
        name = rec[0].rstrip()
        name = name[1:]
        name = name.split()[0]
        seq = rec[1].rstrip()
        seqList = list(seq)
        qualString = rec[3].rstrip()
        lord = ord #local function pointer for speedup 
        qualList = [lord(i) - 33 for i in qualString ]  
        
        for i in range(len(seqList)):
            if qualList[i] <= 2:
               seqList[i] = 'N'
        seqStr = ''.join(seqList)
        outFile.write('>%s\n%s\n' % (name,seqStr))
        n += 1
        if n % 50000 == 0:
            print '\t Did %i records...' % n
    inFile.close()
    outFile.close()
#####################################################################
#def run_exonerate(myData):
#    myData['exonerateOutput'] = myData['outDir'] + 'exonerate.out'
#    
#    cmd = 'exonerate --showcigar yes --showalignment no --showvulgar no --model affine:global -E True %s %s > %s' % (myData['flash_Frag_fasta'], myData['targetFA'],myData['exonerateOutput'])
#    print cmd
#    runCMD(cmd)
#####################################################################
def get_zipcode_noindel(myData):
    myData['extractionRaw'] = myData['outDir'] + 'extraction.raw.txt.gz'
    if os.path.isfile(myData['extractionRaw']) is True:
        print 'Looks like already ran extraction raw'
        return
        
    inFile = gzip.open(myData['flash_Frag'],'r')
    outFile = gzip.open(myData['extractionRaw'],'w')
    n = 0
    print 'Filtering fastq...'
    while True:
        rec = get_4l_record(inFile)
        if rec == '':
            break
        name = rec[0].rstrip()
        name = name[1:]
        name = name.split()[0]
        seq = rec[1].rstrip()
        seqList = list(seq)
        qualString = rec[3].rstrip()
        lord = ord #local function pointer for speedup 
        qualList = [lord(i) - 33 for i in qualString ]  
        
        for i in range(len(seqList)):
            if qualList[i] <= 2:
               seqList[i] = 'N'
        seqStr = ''.join(seqList)
        seqStrRC = revcomp(seqStr)

        maxMatchL = 0
        maxOffsetL = 0
        for offset in range(0,5):
            numMatches = count_matches(seqStr,myData['leftTarget'],offset)
            if numMatches > maxMatchL:
                maxMatchL = numMatches
                maxOffsetL = offset

        maxMatchLRC = 0
        maxOffsetLRC = 0
        for offset in range(0,5):
            numMatches = count_matches(seqStrRC,myData['leftTarget'],offset)
            if numMatches > maxMatchL:
                maxMatchLRC = numMatches
                maxOffsetLRC = offset

        if maxMatchLRC > maxMatchL:
            maxMatchL = maxMatchLRC
            maxOffsetL = maxOffsetLRC
            seqStr = seqStrRC
            
        maxMatchR = 0
        maxOffsetR = 0
        for offset in range(len(seqStr)-len(myData['rightTarget'])-5,len(seqStr)-len(myData['rightTarget'])):
            numMatches = count_matches(seqStr,myData['rightTarget'],offset)
            if numMatches > maxMatchR:
                maxMatchR = numMatches
                maxOffsetR = offset
        
        
        leftPre = seqStr[0:maxOffsetL]
        leftMatch = seqStr[maxOffsetL:maxOffsetL+len(myData['leftTarget'])]                        

        zipCode = seqStr[maxOffsetL+len(myData['leftTarget']):maxOffsetR]

        rightMatch = seqStr[maxOffsetR:maxOffsetR+len(myData['rightTarget'])]
        rightPost = seqStr[maxOffsetR+len(myData['rightTarget']):]
        
        extractionResult = [leftPre,str(maxMatchL),zipCode,str(maxMatchR),rightPost]
        extractionResult = ':'.join(extractionResult) + '\n'
        outFile.write(extractionResult)
        
        

        n += 1
        if n % 50000 == 0:
            print '\t Did %i records...' % n
            
    inFile.close()
    outFile.close()
#####################################################################    
def count_matches(seq1,seq2,seq2offset):  #offset is where in seq1 seq2 starts
    # assumes seq2 is shorter
    len2 = len(seq2)
    s1Part = seq1[seq2offset:seq2offset+len2]    
    numMatches = 0
    numMismatches = 0
    for i in range(len(seq2)):
        if s1Part[i] == seq2[i]:
            numMatches += 1
        else:
            numMismatches += 1
    return numMatches
#####################################################################
def count_extracted_zips(myData):
    myData['zipTable'] = myData['outDir'] + 'extraction.ziptable.gz'

    
    inFile = gzip.open(myData['extractionRaw'],'r')
    
    myData['countsInfo'] = {}
    myData['countsInfo']['failMinLeft'] = 0
    myData['countsInfo']['failMinRight'] = 0    
    myData['countsInfo']['failMinZip'] = 0    
    myData['countsInfo']['failMaxZip'] = 0    
    myData['countsInfo']['zipHasN'] = 0        
    
    myData['countsInfo']['PassZip'] = 0    

    myData['countsInfo']['zipLens'] = {}
    for i in range(myData['minZipLen'],myData['maxZipLen'] +1 ):
        myData['countsInfo']['zipLens'][i] = [0,0.0]   #num zipcodes, percent of reads
        



    zipCodes = {}
    endsPass = {}
    
    for line in inFile:
        line = line.rstrip()
        line = line.split(':')
        leftBC = line[0]
        matchLeft = int(line[1])
        zipCode = line[2]
        matchRight = int(line[3])
        rightBC = line[4]
        
        failReasons = 0
        if matchLeft < myData['minLeftMatch']:
            failReasons += 1
            myData['countsInfo']['failMinLeft'] += 1
        if matchRight < myData['minRightMatch']:
            failReasons += 1
            myData['countsInfo']['failMinRight'] += 1

        if len(zipCode) < myData['minZipLen']:
            failReasons += 1
            myData['countsInfo']['failMinZip'] += 1

        if len(zipCode) > myData['maxZipLen']:
            failReasons += 1
            myData['countsInfo']['failMaxZip'] += 1
            
        if 'N' in zipCode:
            failReasons += 1
            myData['countsInfo']['zipHasN'] += 1
            
            
            
        if failReasons == 0:
            myData['countsInfo']['PassZip'] += 1
            if zipCode not in zipCodes:
                zipCodes[zipCode] = 0
            zipCodes[zipCode] += 1
            
            # check the keys for the ends...
            k = leftBC + '-' + rightBC
            if k not in endsPass:
                endsPass[k] = 0
            endsPass[k] += 1
            
            
    inFile.close()
    
    print 'Number of unique zipcodes:',len(zipCodes)
    print 'Number of unique flanking primer tails for pass set:',len(endsPass)
    
    outFile = gzip.open(myData['zipTable'] ,'w')
    print 'Sorting....',
    zipK = zipCodes.keys()
    zipK.sort(key=lambda k: zipCodes[k],reverse=True)
    print 'Done'
    
    totalDepth = 0
    for k in zipK:
        totalDepth += zipCodes[k]
    print 'Total depth for assigned zipcodes is',totalDepth

    # get top 10 flanks
    primerTails  = endsPass.keys()
    primerTails.sort(key= lambda k: endsPass[k], reverse=True)
    myData['top10primerTails'] = []
    for i in range(10):
        if i >= len(primerTails):
            break
        k = primerTails[i]
        reads = endsPass[k]
        fraction = float(reads)/float(totalDepth)
        myData['top10primerTails'].append([k,reads,fraction])


    for k in zipK:
        outFile.write('%s\t%i\t%.8f\n' % (k,zipCodes[k],float(zipCodes[k])/float(totalDepth)))
        
        myData['countsInfo']['zipLens'][len(k)][0] += 1
        myData['countsInfo']['zipLens'][len(k)][1] += float(zipCodes[k])/float(totalDepth)        
    outFile.close()
#####################################################################
def print_extraction_stats(myData):
    myData['extractStatsFile'] = myData['outDir'] + 'extraction.stats.txt'
    outFile = open(myData['extractStatsFile'],'w')
    
    outFile.write('name\t%s\n' % myData['name'])
    outFile.write('fq1\t%s\n' % myData['fq1'])
    outFile.write('fq2\t%s\n' % myData['fq2'])
    outFile.write('zipregion target\t%s\t%s\n' % (myData['leftTarget'],myData['rightTarget']))
    outFile.write('minLeftMatch\t%i\n' % myData['minLeftMatch'] )
    outFile.write('minRightMatch\t%i\n' % myData['minRightMatch'] )
    outFile.write('zipcode length range\t%i\t%i\n' % (myData['minZipLen'],myData['maxZipLen']))

    outFile.write('\nFLASH read overlap statistics\n')
    outFile.write('total pairs\t%i\n' % myData['flashInfo']['totalpairs'])
    outFile.write('combined pairs\t%i\n' % myData['flashInfo']['combinedpairs'])
    outFile.write('uncombinedpairs pairs\t%i\n' % myData['flashInfo']['uncombinedpairs'])
    outFile.write('fraction combined\t%f\n' % myData['flashInfo']['fraccombined'])
    
    outFile.write('\nZipcode extraction statistics\n')
    klist = ['failMinLeft','failMinRight','zipHasN','failMinZip','failMaxZip']
    for k in klist:
        outFile.write('%s\t%i\n' % (k, myData['countsInfo'][k]))
    outFile.write('Reads with passing zipcode\t%i\n' % (myData['countsInfo']['PassZip']))
    f = float(myData['countsInfo']['PassZip'])/ float(myData['flashInfo']['combinedpairs'])
    outFile.write('Fraction of combined reads with passing zipcode\t%f\n' % f)
    
    outFile.write('\nZipcode size profile\n')
    
    
    klist = myData['countsInfo']['zipLens'].keys()
    klist.sort()
    outFile.write('length\tnumber of zipcodes\tfraction of assigned reads\n')
    for k in klist:
        outFile.write('%i\t%i\t%f\n' % (k,myData['countsInfo']['zipLens'][k][0],myData['countsInfo']['zipLens'][k][1]))

    
    outFile.write('\nPrimer Tail profile\n')
    outFile.write('sequence\tNumber of Reads\tfraction of assigned reads\n')
    for i in myData['top10primerTails']:
        outFile.write('%s\t%i\t%f\n' % (i[0],i[1],i[2]) )
    
    outFile.close()
#####################################################################

