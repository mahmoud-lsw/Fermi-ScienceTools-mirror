/** @file SpectObj.h
  @brief declaration n of SpectObj class

  $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/SpectObj/SpectObj/SpectObj.h,v 1.25 2009/08/06 10:26:08 razzano Exp $

*/
#ifndef SpectObj_H
#define SpectObj_H
#include "TH2D.h"
#include "TH1D.h"
#include "TROOT.h"
#if 0 //THB
#include "TRandom.h"
#endif
#include <vector>
#include <iostream>

/*! 
  \struct photon
  Trivial struct for storing time and energy.
*/
namespace IRB {
   class EblAtten;
}

struct photon
{
  double time;
  double energy;
};

/*! 
  \class SpectObj
  \brief general interface for astronomical spectrum.
  
  SpectObj takes as input a two dimensional histogram containing the photons per square meter as a function of the energy of the time.
  SpectObj compute both the relevant integrals (in energy and time), and both extracts photons to feed the simulators (flux). 

*/
class SpectObj
{
#if 0 //THB
#else
  class UniformRandom {
  public:
    UniformRandom(){}
    double Uniform(double a=0, double b=1);
  };
#endif

 public:
  SpectObj(const TH2D* In_Nv, int type, double z=0.0); //Max
  
  ~SpectObj()
    {
      delete spec;
      delete times;
      delete Probability;
      delete Nv;
      delete m_SpRandGen;
      if (PeriodicSpectrumIsComputed == true ) 
	{
	  delete PeriodicSpectrum;
	  delete PeriodicLightCurve;
	}
      std::cout<<" SpectObj: Generated photons : "<<counts<<" over "<<m_AreaDetector<<" m^2 "<<std::endl;
    }
  void GetUniqueName(void *ptr, std::string & name);
  
  TH1D *Integral_E(double e1, double e2);
  TH1D *Integral_E(int ei1, int ei2); 
  TH1D *Integral_T(double t1, double t2, double e1, double e2);
  TH1D *Integral_T(int ti1, int ti2, int ei1, int ei2);
  TH1D *Integral_T(double t1, double t2, double en = 0.0);
  TH1D *Integral_T(int ti1, int ti2, int ei = 1);
  double Integral_E(TH1* Sp, double e1=0.0, double e2=0.0);
  double Integral_E(TH1* Sp, int ei1, int ei2);
  double Integral_T(TH1* Lc, double t1=0.0, double t2=0.0);
  double Integral_T(TH1* Lc, int ti1, int ti2);
  void ComputeProbability(double enph);
  TH1D *N(TH1D *EN);
  photon GetPhoton(double t0, double enph);
  
  void SetAreaDetector(double AreaDetector=6.0);
  inline void SetFluxFactor(double FluxFactor=1.0) {m_FluxFactor = FluxFactor;};
    
  inline double GetAreaDetector() { return m_AreaDetector; }
  double flux(double time, double enph);
  double interval(double time, double enph);
  double energy(double time, double enph);
  void SaveParameters(double tstart, std::pair<double,double> direction);
  //////////////////////////////////////////////////
  void ScaleAtBATSE(double fluence);
  double GetFluence(double BL=0.0, double BH=0.0);
  double GetPeakFlux(double BL=0.0, double BH=0.0, double AccumulationTime = 0.256);
  double GetT90(double BL=0.0, double BH=0.0);
  void GetGBM();
  //////////////////////////////////////////////////
  TH1D* GetSpectrum(double t=0.0);
  TH1D* GetTimes(double t=0.0);
  TH1D *CloneSpectrum();
  TH1D *CloneTimes();
    
 private:
  int counts;
  double  m_AreaDetector;
  double  m_FluxFactor;

#if 0 //THB
  TRandom *m_SpRandGen;
#else
  UniformRandom * m_SpRandGen;
#endif
  TH2D* Nv;
  int ne,nt;
  int sourceType; //"0=Transient,1=Periodic"
  double emin,emax;

  double m_Tmin,m_Tmax, m_TimeBinWidth;
  double Ptot;
  double m_z;
  double m_meanRate;
  
  TH1D *spec,*times,*Probability,*PeriodicSpectrum,*PeriodicLightCurve;
  photon ph;
  bool ProbabilityIsComputed, PeriodicSpectrumIsComputed;
  IRB::EblAtten * m_tau;

   static bool s_gRandom_seed_set;

};
#endif

