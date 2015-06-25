# $Id: Parfile.py,v 1.1 2008/08/28 16:19:20 elwinter Exp $

# Python module to process HEADAS parameter (.par) files.

#******************************************************************************

# Import external modules.
import os
import re

# Standard modules

# Third-party modules

# Project modules

#******************************************************************************

class Parfile:

    #--------------------------------------------------------------------------

    def __init__(self, filename = None, *args, **kwargs):

        # Find the requested parfile.
        self._path = self._find(filename)

        # Parse the parfile.
        self._parse()

    #--------------------------------------------------------------------------

    def GetParameterValue(self, name):
        return self._parameters[name]['value']

    #--------------------------------------------------------------------------

    def _find(self, filename):

        # Search the directories listed in the PFILES environment
        # variable for the requested parfile. Return the full path if
        # found. Raise an IOError exception if not found.
        parfile = ''
        for path in os.environ['PFILES'].split(';'):
            parfile = path + os.sep + filename
            if os.path.exists(parfile):
                return parfile
        raise IOError

    #--------------------------------------------------------------------------

    def _parse(self):
        self._parameters = {}
        for line in open(self._path):

            # Skip blank lines.
            if re.match('^\s*$', line):
                continue

            # Skip comment lines.
            if re.match('^\s*#', line):
                continue

            # Strip the trailing newline.
            line = line.rstrip()

            # Split the line on commas.
            fields = line.split(',')

            # Save the fields for this record.
            name = fields[0]
            type = fields[1]
            field2 = fields[2]
            value = fields[3].strip('"')
            minimum = fields[4]
            maximum = fields[5]
            comment = fields[6].strip('"')
            parameter = { 'type' : type,
                          'field2' : field2,
                          'value': value,
                          'minimum': minimum,
                          'maximum': maximum,
                          'comment': comment }
            self._parameters[name] = parameter

#******************************************************************************

# Self-test code.

# If this code generates any output, an error has been found.

if __name__ == '__main__':
    parfile = Parfile('gtobspsf.par')
    print parfile._path
