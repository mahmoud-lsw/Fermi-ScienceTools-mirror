"""
@brief Pass 4 event classes based on Toby's straw man proposal.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/fitsGen/python/Pass4_Classifier.py,v 1.3 2007/04/08 19:49:03 jchiang Exp $
#
from EventClassifier import EventClassifier

meritVariables = """
EvtRun    EvtEnergyCorr 
McEnergy  McXDir  McYDir  McZDir   
McXDirErr   McYDirErr  McZDirErr   
McTkr1DirErr  McDirErr  
GltWord   FilterStatus_HI 
Tkr1FirstLayer  
CTBCORE  CTBSummedCTBGAM  CTBBestEnergyProb
""".split()

#
# Event class cuts for standard/front, standard/back
#

eventClassCuts = ['&&'.join(('(CTBBestEnergyProb>0.1)', '(CTBCORE>0.1)',
                             '(CTBSummedCTBGAM>0.5)', '(Tkr1FirstLayer>5.5)')),
                  '&&'.join(('(CTBBestEnergyProb>0.1)', '(CTBCORE>0.1)',
                             '(CTBSummedCTBGAM>0.5)', '(Tkr1FirstLayer<5.5)'))]

eventClassifier = EventClassifier(eventClassCuts)
