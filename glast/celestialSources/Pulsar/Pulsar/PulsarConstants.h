/**
 * @file PulsarConstant.h
 * @brief Collection of constants useful for PulsarSpectrum simulator
 * @ author Massimiliano Razzano (massimiliano.razzano@pi.infn.it
 * @ author Nicola Omodei (nicola.omodei@pi.infn.it
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/Pulsar/Pulsar/PulsarConstants.h,v 1.36 2008/07/22 16:46:49 razzano Exp $
 */
#ifndef PulsarCONSTANT_HH

#define  PulsarCONSTANT_HH 



#include <vector> 

#include <cmath>

#include <string>

#include "TRandom.h"



/*! 

 * \class PulsarConstants

 * \brief Class instantiated to access general parameters and constants.

 *  

 * The namespace cst contains all the constant needed to the simulation.

 *

 *

 * \author Nicola Omodei        nicola.omodei@pi.infn.it 

 * \author Massimiliano Razzano massimiliano.razzano@pi.infn.it

 */





//! Namespace containing general parameters and constants

namespace cst

{

  //! Lower Extremum of the normalization band,espressed in keV (100MeV) 

  const double EnNormMin = 1e5;  



  //! Upper Extremum of the normalization band,espressed in keV (30GeV) 

  const double EnNormMax = 3e7;  

  

  //! Number of energy bins used in the TH2D histogram

  const    int Ebin =  100 ; 



  //! Number of time bins used in the TH2D histogram

  const int Tbin =  200; 


  //! Conversion erg-->MeV

  const double erg2meV   = 624151.0;



  //! Spedd of light (km/s)

  const double clight = 299792.45;



  //! Lower energy of GBM band, expressed in keV (10keV)

  const double GBM1=10.0;                     

  

  //! Upper energy of GBM band, expressed in keV (25MeV)

  const double GBM2=2.5e4;                   



  //! Lower energy of LAT band, expressed in keV (30MeV)

  const double LAT1=3.0e4;                     



  //! Upper energy of LAT band, expressed in keV (300GeV) 

  const double LAT2=3.0e8;                     



  //! Lower energy of EGRET band, expressed in keV (30MeV)

  const double EGRET1=3.0e4;                 



  //! Mid energy of EGRET band, expressed in keV (100MeV)

  const double EGRET2=1.0e5;                  



  //! Upper energy of EGRET band, expressed in keV (30GeV)

  const double EGRET3=3.0e7;                  



  //! Start Mission Date, expressed in MJD (corresponding to TT January 1st, 2007, at 00:00:00)

  const double StartMissionDateMJD = 51910.0 + 64.184/86400.;//54101.0; 

  

  //! Tolerance for the ephemerides decorrection (in s.)

  const double ephemCorrTol = 1e-6;



  //! Tolerance for the barycentri decorrection (in us.)

  const double baryCorrTol = 5e-6;



  //! Tolerance for the inverse binary demodulation (in s.)

  const double InverseDemodTol = 5e-6;



  //! Tolerance for the iterative binary demodulation

  const double DemodTol = 1e-8;



  //! Difference between JD and MJD

  const double JDminusMJD = 2400000.5; 



  //! Seconds in one day

  const int SecsOneDay = 86400; 



  //! Degree to Radiants conversion

  const double DegToRad = M_PI/180.;



}



#endif

