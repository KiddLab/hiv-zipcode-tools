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

#####################################################################
# setup paths to default programs to use....
def set_default_prog_paths(myData):
    print 'Checking for required programs\n'

    
    print 'Checking flash...'
    if which('flash') is None:
        print 'flash not found in path! please fix (module load?)'
        sys.exit()

    print 'Checking exonerate...'
    if which('exonerate') is None:
        print 'exonerate not found in path! please fix (module load?)'
        sys.exit()
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
    
    if os.path.isfile(myData['flash_1']) is False:
        cmd = 'flash'  +' -r 75 -f 70 -s 2 --cap-mismatch-quals -z -t 1 ' + ' %s %s -d %s | tee %s' % (myData['fq1'],myData['fq2'],myData['flashDir'],myData['flash_stats'])
        print cmd
        runCMD(cmd)
    else:
        print 'looks like flash already ran!'
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
    myData['flash_Frag_fasta'] = myData['flashDir'] + 'out.extendedFrags.fasta.gz'
    if os.path.isfile(myData['flash_Frag_fasta']) is True:
        print 'Looks like already ran filter fasta'
        return
        
    inFile = gzip.open(myData['flash_Frag'],'r')
    outFile = gzip.open(myData['flash_Frag_fasta'],'w')
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



#    record['qualString'] = record['qualString'].rstrip()    
#    lord = ord #local function pointer for speedup    
#    record['qualList'] = [lord(i) - qoffSet for i in record['qualString'] ]  
#    record['qual33Str'] = ''.join([chr(i+33) for i in record['qualList']])




