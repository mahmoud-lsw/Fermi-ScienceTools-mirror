"""
This module uses GtApp to wraps the Science Tools as python objects.
"""
#
# $Header: /glast/ScienceTools/glast/sane/python/gt_apps.py,v 1.1.1.3 2012/01/17 17:18:11 elwinter Exp $
#
from GtApp import GtApp

#
# Likelihood applications
#
expCube = GtApp('gtltcube', 'Likelihood')
addCubes = GtApp('gtltsum', 'Likelihood')
expMap = GtApp('gtexpmap', 'Likelihood')
diffResps = GtApp('gtdiffrsp', 'Likelihood')
like = GtApp('gtlike', 'Likelihood')
TsMap = GtApp('gttsmap', 'Likelihood')

counts_map = GtApp('gtbin', 'evtbin')
srcMaps = GtApp('gtsrcmaps', 'Likelihood')
model_map = GtApp('gtmodel', 'Likelihood')
gtexpcube2 = GtApp('gtexpcube2', 'Likelihood')

#
# observationSim
#
obsSim = GtApp('gtobssim', 'observationSim')
#orbSim = GtApp('gtorbsim', 'observationSim')

#
# dataSubselector
#
filter = GtApp('gtselect', 'dataSubselector')
maketime = GtApp('gtmktime', 'dataSubselector')

##
## Pulsar-related
##
#pulsePhase = GtApp('gtpphase', 'pulsePhase')
#psearch = GtApp('gtpsearch', 'periodSearch')

#
# others
#
evtbin = GtApp('gtbin', 'evtbin')
rspgen = GtApp('gtrspgen', 'rspgen')



