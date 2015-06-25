#include "GRBobs/GRBobsPulse.h"
#include "GRBobs/GRBobsConstants.h"
#include <iostream>
#include <cmath>

#define DEBUG 0

using std::fabs; using std::pow;

GRBobsPulse::GRBobsPulse(){;}

GRBobsPulse::GRBobsPulse(double peakTime, GRBobsParameters *m_params)
{  
  m_peakTime     = peakTime;
  m_riseTime     = m_params->GetRiseTime();
  m_decayTime    = m_params->GetDecayTime();
  m_Peakedness   = m_params->GetPeakedness();
  
  m_Intensity    = m_params->GetPulseHeight();
    
  m_start        = m_peakTime - m_riseTime  * pow(log(100.0),1.0/m_Peakedness);
  m_end          = m_peakTime + m_decayTime * pow(log(100.0),1.0/m_Peakedness);
  
  m_Epeak        = m_params->GetEpeak();
  m_LowEnergy    = m_params->GetLowEnergy();
  m_HighEnergy   = m_params->GetHighEnergy();
  m_duration     = m_end-m_start;
  
  if(DEBUG)  Print();
  
}

void GRBobsPulse::Print()
{
  std::cout<<"Pulse: start= "<<m_start
	   <<" peak= "<<m_peakTime
	   <<" end= "<<m_end
	   <<" rise t= "<<m_riseTime
	   <<" Decay t= "<<m_decayTime
	   <<" Dur.= "<<m_duration
	   <<" Intensity= "<<m_Intensity
	   <<" Peakedness= "<<m_Peakedness
	   <<" Peak E= "<<m_Epeak
	   <<" Low idx= "<<m_LowEnergy-ObsCst::We
	   <<" High idx= "<<m_HighEnergy-ObsCst::We
	   <<std::endl;
}

double GRBobsPulse::PulseShape(double t, double e)
{
  //  double tp = m_ts + m_tp;
  double rt = m_riseTime  * pow(e/ObsCst::E0,-ObsCst::We);
  double dt = m_decayTime * pow(e/ObsCst::E0,-ObsCst::We);
  double deltaTP = ObsCst::deltaTPeak * (m_riseTime - rt) * pow(log(100.),1.0/m_Peakedness);
  
  double pt = m_peakTime - deltaTP;
  double pulse=0;
  if (t<pt)
    {
      pulse = exp(-pow(fabs(t-pt)/rt,m_Peakedness));
    }
  else
    {
      pulse = exp(-pow(fabs(t-pt)/dt,m_Peakedness));
    }
  
  //////////////////////////////////////////////////
  double a  = m_LowEnergy;
  double b  = m_HighEnergy;
  
  double Ep = m_Epeak;
  double E0 = Ep/(2.+a-ObsCst::We); 
  double Ec = (a-b)*E0;
  double C  = pow(Ec/100.,a-b)*exp(b-a);
  
  double bandf;
  if(e < Ec) 
    bandf = pow(e/100.,a) * exp(-e/E0);
  else
    bandf= C * pow(e/100.,b); // ph cm^(-2) s^(-1) keV^(-1)

  if(DEBUG)
    {
      std::cout<<" Parameters used for the band functions:"<<std::endl;
      std::cout<<" a = "<<a<<" b= "<<b<<" Ep= "<<Ep<<" E0= "<<E0<<" Ec= "<<Ec<<std::endl;
    }


  return m_Intensity * bandf * pulse;
}



