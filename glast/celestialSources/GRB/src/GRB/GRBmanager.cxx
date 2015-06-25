#include "GRBmanager.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 
#include "astro/GPS.h"
#include "astro/SkyDir.h"
#include "astro/PointingTransform.h"
#include "astro/JulianDate.h"

#include "CLHEP/Vector/ThreeVector.h"

ISpectrumFactory &GRBmanagerFactory() 
{
  static SpectrumFactory<GRBmanager> myFactory;
  return myFactory;
}


GRBmanager::GRBmanager(const std::string& params)
  : m_params(params)
{
  m_Nbursts=0;
  paramFile = facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("GRB"),"GRBParam.txt");

  facilities::Util::expandEnvVar(&paramFile);
  
  m_startTime       = parseParamList(params,0)+Spectrum::startTime();//m_startTime   = TMath::Max(0.,parseParamList(params,0));
  m_timeToWait  = TMath::Max(0.,parseParamList(params,1));
  m_enph        = TMath::Max(0.,parseParamList(params,2))*1000.0; //keV

  m_par = new Parameters();
  //////////////////////////////////////////////////
  GenerateGRB();  
  //////////////////////////////////////////////////
  if(m_enph<=0) m_enph=cst::enph;
}

GRBmanager::~GRBmanager() 
{  
  //  std::cout<<"~GRBmanager() "<<std::endl;
  delete m_par;
  delete m_GRB;
  delete m_spectrum;
}

//return flux, given a time
double GRBmanager::flux(double time) const
{


  if(time <= m_startTime || (time > m_endTime)) 
    return 0.0;
  return m_spectrum->flux(time-m_startTime,m_enph);
}

TString GRBmanager::GetGRBname(double time)
{
  ////DATE AND GRBNAME
  /// GRB are named as customary GRBYYMMDDXXXX
  //  The mission start and launch date are retrieved, 
  //  and the current burst's start time is added to them to form a Julian Date.
  //  Note: the argument of the JulianDate constructor is in day.
  
  astro::JulianDate JD(astro::JulianDate::missionStart() +
		       (Spectrum::startTime()+m_startTime)/astro::JulianDate::secondsPerDay);

  int An,Me,Gio;
  
  m_Rest=((int) m_startTime % 86400) + m_startTime - floor(m_startTime);
  m_Frac=(m_Rest/86400.);
  int FracI=(int)(m_Frac*1000.0);
  double utc;
  JD.getGregorianDate(An,Me,Gio,utc);

  TString GRBname="";

  An-=2000;
  
  if(An<10) 
    {
      GRBname+="0";
      GRBname+=An;
    }
  else
    {
      GRBname+=An;
    }
  if(Me<10) 
    {
      GRBname+="0";
      GRBname+=Me;
    }
  else
    {
      GRBname+=Me;
    }
  if(Gio<10) 
    {
      GRBname+="0";
      GRBname+=Gio;
    }
  else
    {
      GRBname+=Gio;
    }


  if (FracI<10)
    {
      GRBname+="00";
      GRBname+=FracI;
    }
  else if(FracI<100) 
    {
      GRBname+="0";
      GRBname+=FracI;
    }
  else 
    GRBname+=FracI;
    
  //std::cout<<"GENERATE GRB ("<<GRBname<<")"<<std::endl;
  
  //..................................................  //
  return GRBname;
}


void GRBmanager::GenerateGRB()
{
  //  Field Of View for generating bursts degrees above the XY plane.
  const double FOV = -100;
  if(FOV>-90) std::cout<<"WARNING!! GRB: FOV = "<<FOV<<std::endl;

  using astro::GPS;
  //////////////////////////////////////////////////
  m_Nbursts++;
  m_par->ComputeParametersFromFile(paramFile,m_Nbursts);
  m_GRB      = new  GRBSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->Fireball(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  m_nextBurst = m_endTime   + m_par->rnd->Exp(m_timeToWait);
  m_fluence   = m_GRB->GetFluence();
  m_GRBnumber = m_GRB->GetGRBNumber();

  m_theta = -10000000.0;
  while(m_theta<FOV)
    {
      m_par->SetGalDir(-200,-200); //this generates random direction in the sky
      m_GalDir = m_par->GetGalDir();      
      m_l = m_GalDir.first;
      m_b = m_GalDir.second;
      
      astro::SkyDir sky(m_l,m_b,astro::SkyDir::GALACTIC);
      CLHEP::Hep3Vector skydir=sky.dir();
      CLHEP::HepRotation   rottoglast = GPS::instance()->transformToGlast(m_startTime,GPS::CELESTIAL);
      CLHEP::Hep3Vector scdir = rottoglast * skydir;
      m_ra    = sky.ra();
      m_dec   = sky.dec();
      m_theta = 90. - scdir.theta()*180.0/M_PI; // theta=0 -> XY plane, theta=90 -> Z
      m_phi   = scdir.phi()*180.0/M_PI;
    }
  TString GRBname = GetGRBname(m_startTime);
  std::ofstream os;
  if(m_Nbursts==1) 
    {
      os.open("grb_generated.txt",std::ios::out);
      os<<"m_GRBnumber   GRBname   m_startTime   m_endTime   m_ra   m_dec   m_l    m_b   m_theta   m_phi   m_fluence  m_ic"<<std::endl;
    }
  else 
    os.open("grb_generated.txt",std::ios::app);
  
  os<<m_GRBnumber<<" "<<GRBname<<" "<<m_startTime<<" "<<m_endTime;
  os<<" "<<m_ra<<" "<<m_dec<<" "<<m_l<<" "<<m_b<<" "<<m_theta<<" "<<m_phi;
  os<<" "<<m_fluence<<" "<<m_par->GetInverseCompton()<<std::endl;

  os.close();

  std::cout<<"Physical Model GRB"<<GRBname<<" t start "<<m_startTime<<", tend "<<m_endTime
	   <<" l,b = "<<m_l<<", "<<m_b<<" elevation,phi(deg) = "<<m_theta<<", "<<m_phi<<" Fluence = "<<m_fluence<<std::endl;
  
  if(m_par->GenerateGBM())
    {
      m_GRB->SaveGBMDefinition(GRBname,m_ra,m_dec,m_theta,m_phi,m_startTime);
      m_GRB->GetGBMFlux(GRBname);
    }
}

double GRBmanager::interval(double time)
{ 
  double inte;  
  if(time <= m_startTime) 
    inte = m_startTime - time + m_spectrum->interval(0.0,m_enph);
  else if (time<m_endTime)
    inte = m_spectrum->interval(time - m_startTime,m_enph);
  else 
    {
      delete m_GRB;
      delete m_spectrum;
      m_startTime = m_nextBurst;
      GenerateGRB(); 
      inte = m_startTime-time + m_spectrum->interval(0.0,m_enph);
    }
  return TMath::Min(inte,m_nextBurst-time);
}

double GRBmanager::energy(double time)
{
  return m_spectrum->energy(time-m_startTime,m_enph)*1.0e-3; //MeV
}

double GRBmanager::parseParamList(std::string input, int index)
{
  std::vector<double> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    double f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  if(index>=(int) output.size()) return 0.0;
  return output[index];
}


