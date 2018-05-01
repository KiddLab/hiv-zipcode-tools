# zipcodetools.py
# set of functions for working with revised HIV zipcodes
# python 2.X code
import sys
import genutils
import os
import gzip

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
#####################################################################
# setup paths to default programs to use....
def set_default_prog_paths(myData):
    print 'Checking flash...'
    if which('flash') is None:
        print 'flash not found in path! please fix (module add?)'
        sys.exit()

    print 'Checking exonerate...'
    if which('exonerate') is None:
        print 'exonerate not found in path! please fix (module add?)'
        sys.exit()
#####################################################################



