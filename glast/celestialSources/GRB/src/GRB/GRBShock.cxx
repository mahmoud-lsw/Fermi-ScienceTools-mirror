//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//
#include <cmath>
#include <iostream>

#include "GRBShock.h"
#include "GRBShell.h"

#define DEBUG 0

#define PulseShape 0 

#define Spectrum 0 

using namespace cst;
//////////////////////////////////////////////////
// BS --> FS ->
//BS = Fast Shell
//FS = Slow Shell
GRBShock::GRBShock(GRBShell *FS, GRBShell *BS, double tshock, double p)
{
  m_p  =    p;
  m_IC =  1.0;
    
  double g1 = FS->GetGamma();
  double g2 = BS->GetGamma();
  
  //  double e1 = FS->GetEnergy();
  double r2 = BS->GetRadius();

  double dr1 = FS->GetThickness();
  double dr2 = BS->GetThickness();

  double m1 = FS->GetMass();
  double m2 = BS->GetMass();
  
  tsh = tshock - r2/c;
  
  double M,b12,b1,b2;
  
  // Initial energy of the final shell: (obs)
  ei  = (m1*g1+m2*g2) * c2; //erg
  
  // mass of the final shell:
  mf  = m1+m2; //g
  
  //  std::cout<<"--------------------My Model--------------------"<<std::endl;
  M   = sqrt(m1*m1+m2*m2+2.0*m1*m2*(g1*g2-sqrt(g1*g1-1.0)*sqrt(g2*g2-1.0)))-(mf);
  // Final: (obs) 
  gf  = (m1*g1 + m2*g2)/(mf+M);
  
  //gf  = sqrt((m1*g1 + m2*g2)/(m1/g1 + m2/g2));
  ef  = mf * gf * c2;       //erg
  // Internal : (obs) 
  eint_o = ei - ef;         //erg
  eint_c = (eint_o)/gf;     //erg
  eff = eint_o/ei;
  //////////////////////////////////////////////////
  // Radius of the shell:
  rf  = r2;//cm

  // Thickness of the final shell
  drf = (dr2 + dr1)/2.;//cm
  MS = new GRBShell(gf,rf,drf,ef);
  
  b1  = FS->GetBeta();
  b2  = BS->GetBeta();
  // relative velocity:
  b12 = (b1-b2)/(1.0-b1*b2);
  double g12 = 1.0/sqrt(1.0 - b12*b12);
  double n1  = FS->GetComPartDens(); //1/cm^3
 
  //  double n2  = BS->GetComPartDens(); //1/cm^3
  //  double x = n1/n2;
  // Lorentz factor of the shock in the comoving frame of the unshocked fluid
  //  gsh =  TMath::Max(1.0,pow(x,1./4.)*sqrt(g12));
  gsh =  sqrt((pow(g12,2.0)+1.0)/2.);
  double bsh=TMath::Max(1.0e-6,sqrt(1.0-1.0/(gsh*gsh)));
  // density of the shocked fluid
  nsh  = 4.0 * n1 * pow(gsh,2.0);//1/cm^3
  esh  = nsh * mpc2 / erg2meV; //erg
  // Magnetic field in gauss
  B    = sqrt(8.0 * pi * ab * esh); 
  //  B    = 0.39 * sqrt(ab * nsh);
  
  // synchrotron cooling time for an electron with gamma=1
  ts0  = (5.1e8)/pow(B,2.0)/(2.0*gf); //oss
  //   tcom = 2 gf toss 	   
  // Cooling & time scales in observer frame:
  // a) angular time scale
  tar   =   rf/c;
  ta      =   tar/(2.* gf * gf);  // oss
  // a) Shock crossing time scale:
  tc    =  TMath::Min(ta,drf/(2.0*bsh*c));// oss (2.* c * gf);// Kobayashi, Sari, Piran 97;
  // c) cooling time scale:
  //characteistic gamma factors:
  // a) Minimum:
  gem  = TMath::Max(1.0, (m_p-2.0)/(m_p-1.0)* mpc2/mec2 * ae/csi);
  // b) Cooling :
  gec  = TMath::Max(1.0,ts0/tc);
  // c) Maximum:
  geM  = TMath::Max(1.0,3.77e7/sqrt(B)); 
  
  if(DEBUG)
    {
      if(gem<gec)
	std::cout<<" SLOW COOLING (gem = "<<gem<<" gec = "<<gec<<")"<<std::endl;
      else 
	std::cout<<" FAST COOLING (gem = "<<gem<<" gec = "<<gec<<")"<<std::endl;
      Print();
    }
}

void GRBShock::SetTime(double time)
{
  tsh=time;
}

Double_t GRBShock::Peak(Double_t time,Double_t energy)
{
  using std::fabs; using std::pow;
  Double_t to  = time-tsh;   // oss
  if(PulseShape==0)
    {
      //////////////////////////////////////////////////
      if(to<=0) return 0.0;
      double mu      = TMath::Max(-1.0,1.0 - to/tar);
      double beta    = sqrt(1.0-1.0/pow(gf,2.0));
      double doppler = gf*(1.0-beta*mu);
      double sintheta = sqrt(1.0-mu*mu);
      double peak = sintheta/pow(doppler,4.0);
      return eff*peak;
    }
  //////////////////////////////////////////////////
  double ts    = ts0/GammaElectron(energy); //oss
  
  double width = sqrt(ta*ta + tc*tc + ts*ts);//GammaElectron(energy); //oss;
  //width*=pow(energy,-0.4);
  double ratio = 0.4;
  double decay = width /(1.+ratio);
  double rise  = ratio*decay; 
      
  if(PulseShape==1)
    {
      if(to <= tc)
	return  eff*exp(-pow(fabs(to-tc)/rise,1.0));
      return  eff*exp(-pow(fabs(to-tc)/decay,1.0));
    }
  
  if(to<0) return 0.0;
  Double_t F1 = 1./pow(1.+    to/width,2.);
  Double_t F2 = 1./pow(1.+(to - tc)/width,2.);
  if(to <= tc)
    return  eff*(1.- F1);
  return  eff*(F2 - F1);
}
//////////////////////////////////////////////////
Double_t GRBShock::SynSpectrum(Double_t time, Double_t energy)
{
  double to = time-tsh;
  double ec,em,eM,Fv;

  //////////////////////////////////////////////////
  
  double mu      = TMath::Max(-1.0,1.0 - to/tar);
  double beta    = sqrt(1.0-1.0/pow(gf,2.0));
  double doppler = gf*(1.0-beta*mu);
  energy*=doppler; //como energy
  em = EsynCom(gem);
  ec = EsynCom(gec);
  eM = EsynCom(geM);
  
  if(ec <= em) //FAST COOLOING REGIME
    {
      if(energy<=ec)
	Fv = pow(energy/ec,1./3.);
      else if(energy<em)
	Fv = pow(energy/ec,-1./2.);
      else
	Fv = pow(em/ec,-1./2.)*pow(energy/em,-m_p/2.);
    }
  else  // SLOW COOLING REGIME:
    {
      if(energy<=em)
	Fv = pow(energy/em,1./3.);
      else if(energy<ec)
	Fv = pow(energy/em,-(m_p-1.)/2.);
      else
	Fv = pow(ec/em,-(m_p-1.)/2.)*pow(energy/ec,-m_p/2.);
    }
  
  return Fv * exp(-energy/eM) * Peak(time,energy);  // Fv in [keV/(cm� s keV)]
}

//////////////////////////////////////////////////
Double_t GRBShock::ICSpectrum(Double_t time, Double_t energy)
{ 
  return SynSpectrum(time,energy/(gem*gem)) / (gem*gem) * exp(-energy/(gem*gf*mec2*1.0e3));
}

double GRBShock::ComputeFlux(double time, double energy)
{
  if(time<=tsh) return 0.0;
  //this has to return Nv in [ph/(cm� s keV)]
  return (SynSpectrum(time,energy) + m_IC * ICSpectrum(time,energy))/energy;
}

void GRBShock::Print()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Final Gamma :"<<gf<<" Mass :"<<mf<<" Radius = "<< rf<<" Whidth: "<<drf<<std::endl;
  std::cout<<"t shock = "<<tsh<<" ta = "<<ta<<" tc = "<<tc<<" ts = "<<ts0<<" Duration = "<<GetDuration()<<std::endl;
  std::cout<<"Efficiency: "<<eff<<" Initial E = "<<ei<<" Final E = "<<ef<<std::endl;
  std::cout <<" Int E (obs): "<<eint_o<<" (com): "<<eint_c<<" B Field = "<<B<<std::endl;
  std::cout<<" Nsh = "<<nsh<<" Gsh = "<<gsh<<" p "<<m_p<<std::endl;
  std::cout<<"Gmin = "<<gem<<" Gcool = "<<gec<<" Gmax = "<<geM<<std::endl;
  std::cout<<"Es0 = "<<EsynObs(1)<<" Esm = "<<EsynObs(gem)<<" Esc = "<<EsynObs(gec)<<" EsM = "<<EsynObs(geM)<<" IC/Syn = "<<m_IC<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
}
