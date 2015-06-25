#!/usr/bin/env python
"""
Search the PFILES environment variable for the first instance of the
desired .par file.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Header: /glast/ScienceTools/glast/sane/python/pfilesPath.py,v 1.3 2011/03/22 00:32:09 elwinter Exp $
#

import os, sys, re
#from facilities import py_facilities
#os_environ = py_facilities.commonUtilities_getEnvironment

ParFileError = 'ParFileError'
def pfilesPath(parfile=None):
    if os.environ['FERMI_DIR']:
        path = os.sep.join([os.environ['FERMI_DIR'], 'syspfiles'])
        return path
    try:
        if os.name == 'posix':
            pattern = re.compile(";*:*")
        else:
            pattern = re.compile("#*;*")
        paths = re.split(pattern, os.environ['PFILES'])
#        paths = re.split(pattern, os_environ('PFILES'))
        paths = [path for path in paths if path != '']
    except KeyError:
        print "Your PFILES environment variable is not set."
        raise KeyError
    if parfile is None:
        return paths[0]
    for path in paths:
        try:
            if parfile in os.listdir(path):
                return path
        except OSError:
            pass
    raise ParFileError, ".par file " + parfile + " not found."

if __name__ == "__main__":
    if len(sys.argv) == 2:
        parfile = sys.argv[1]
    else:
        parfile = 'likelihood.par'
    print pfilesPath(parfile)
