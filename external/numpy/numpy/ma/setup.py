#!/usr/bin/env python
from __future__ import division, print_function

__author__ = "Pierre GF Gerard-Marchant ($Author: jasercio $)"
__version__ = '1.0'
__revision__ = "$Revision: 1.3.6.1 $"
__date__     = '$Date: 2015/02/27 05:02:49 $'

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('ma', parent_package, top_path)
    config.add_data_dir('tests')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict()
    setup(**config)
