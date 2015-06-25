/////////////////////////////////////////////////
// File PulsarSpectrum.cxx
// Implementation of PulsarSpectrum class
//////////////////////////////////////////////////

#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include "facilities/commonUtilities.h"
#include "facilities/Util.h"
#include <cmath>
#include <cstdlib>
#include <vector>

#define DEBUG 0

using namespace cst;

//======================================================================
// copied from atErrors.h
//----------------------------------------------------------------------
/* Error Code of atFunctions. */
#define NOT_CONVERGED 10        /*equation was not solved*/
//----------------------------------------------------------------------
// copied from atKepler.c (& modified)
//----------------------------------------------------------------------
//#include "atFunctions.h"
//#include "atError.h"
//#include <math.h>

/*
 * solve Kepler equation (KEPLER)  g + e sin E = E
 */
int atKepler(
        double g,        /* input: mean anomaly */
        double eccent,        /* input: eccentricity */
        double *e)        /* output: eccentric anomaly */
{
    static double eps = 1e-7;
    static int imax = 50;

    int i;
    static double error, deltae, d__1;

    *e = g;
    if (g == 0.) return 0;

    for (i=0; i<imax; i++) {
        deltae = (g - *e + eccent * std::sin(*e)) / (1. - eccent * std::cos(*e));
        *e += deltae;
        error = (d__1 = deltae / *e, std::fabs(d__1));
        if (error < eps) return 0;
    }
    return NOT_CONVERGED;
}

/////////////////////////////////////////////////
ISpectrumFactory &PulsarSpectrumFactory() 
 {
   static SpectrumFactory<PulsarSpectrum> myFactory;
   return myFactory;
 }

/////////////////////////////////////////////////
/*!
 * \param params string containing the XML parameters
 *
 * This method takes from the XML file the name of the pulsar to be simulated, the model
 * to be used and the parameters specific for the choosen model. Then extracts from the
 * TXT DataList file the characteristics of the choosen pulsar (e.g. period, flux, period derivatives
 * ,ephemerides, flux, etc...)
 * The names of DataList files are defined in the file <i>PulsarDataList.txt</i><br>
 * The parameters are:<br>
 * <ul>
 *  <li> Pulsar name;
 *  <li> Flux, expressed in ph/cm2/s (E>100MeV);
 *  <li> Ephemerides type (P=period and derivatives;F=frequency and derivatives)
 *  <li> Period,pdot,p2dot or f0,f1,f2;
 *  <li> Start of ephemerides validity range (in MJD);
 *  <li> Reference Epoch, t0 (in MJD);
 *  <li> End of ephemerides validity range (in MJD);
 * </ul> 
 * The parameters taken from XML file are:<br>
 * <ul>
 *  <li> Pulsar name;
 *  <li> Position (RA,dec) in degrees;
 *  <li> Minimum and maximum energy of extracted photons (in keV);
 *  <li> Model choosen (default = 1, phenomenological);
 *  <li> Random seed;
 *  <li> Model-dependend parameters;
* </ul>
 * Then it calls the PulsarSim class in order to have the TH2D histogram.
 * For more informations and for a brief tutorial about the use of PulsarSpectrum you can look at:
 * <br>
 * <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpTutor/PulsarSpTutor.htm </a>
 */


PulsarSpectrum::PulsarSpectrum(const std::string& params)
  : m_solSys(astro::SolarSystem::EARTH), m_params(params) 
{
  m_ff=false;
  // Reset all variables;
  m_RA = 0.0;
  m_dec = 0.0;
  m_l = 0.0;
  m_b = 0.0;
  m_flux = 0.0;
  m_ephemType = "";
  m_period = 0.0;
  m_pdot = 0.0;
  m_p2dot = 0.0;
  m_f0 = 0.0;
  m_f1 = 0.0;
  m_f2 = 0.0;
  m_f0NoNoise = 0.0;
  m_f1NoNoise = 0.0;
  m_f2NoNoise = 0.0;
  m_N0 = 0.0;
  m_t0Init = 0.0;
  m_t0 = 0.0;
  m_t0End = 0.0;  
  m_phi0 = 0.0;
  m_model = 0;
  m_seed = 0;
  m_TimingNoiseModel = 0;
  m_TimingNoiseActivity = 0.;
  m_TimingNoiseRMS = 0.;
  m_TimingNoiseMeanRate=0.;
  m_TimingNoiseTimeNextEvent=0.;
  m_BinaryFlag = 0;
  m_Porb = 0;
  m_Porb_dot = 0.;
  m_asini = 0;
  m_xdot = 0.;
  m_ecc = 0;
  m_ecc_dot = 0.;
  m_omega = 0;
  m_omega_dot = 0.;
  m_gamma = 0.;
  m_shapiro_r = 0.;
  m_shapiro_s = 0.;
  m_t0PeriastrMJD = 0;
  m_t0AscNodeMJD = 0;
  m_PPN =0.;
  m_FT2_startMET = -1.;
  m_FT2_stopMET =-1.; 
  m_Sim_startMET = -1.;
  m_Sim_stopMET =-1.; 

  m_UseFT2 = 0;

  //Read from XML file
  m_PSRname    = parseParamList(params,0).c_str();            // Pulsar name
  m_RA = std::atof(parseParamList(params,1).c_str());         // Pulsar Right Ascension
  m_dec = std::atof(parseParamList(params,2).c_str());         // Pulsar Declination
  m_enphmin    = std::atof(parseParamList(params,3).c_str()); // minimum energy of extracted photons
  m_enphmax    = std::atof(parseParamList(params,4).c_str()); // minimum energy of extracted photons
  m_model      = std::atoi(parseParamList(params,5).c_str()); // choosen model
  m_seed       = std::atoi(parseParamList(params,6).c_str()); //Model Parameters: Random seed

  //Setting random seed
  m_PSpectrumRandom = new TRandom;
  m_PSpectrumRandom->SetSeed(m_seed);

 
  if (m_model == 1) //Phenomenological model
    {
      m_ppar0   = std::atoi(parseParamList(params,7).c_str()); //Model Parameters: Number of peaks  
      m_ppar1 = std::atof(parseParamList(params,8).c_str());   // model parameters
      m_ppar2 = std::atof(parseParamList(params,9).c_str());
      m_ppar3 = std::atof(parseParamList(params,10).c_str());
      m_ppar4 = std::atof(parseParamList(params,11).c_str());
    }
  else 
    if (m_model == 2) //PulsarShape model
      {
	m_ppar0   = std::atoi(parseParamList(params,7).c_str());  //Model Parameters: Use normalization?
	m_PSRShapeName = parseParamList(params,8).c_str();        // model parameters
      }

  //look for pulsar data directory, following, in the order $PULSARDATA, $SKYMODEL_DIR/pulsars or $PULSARROOT/data
  //  std::string m_pulsardata_dir;

  if (::getenv("PULSAR_OUTPUT_LEVEL"))
    {
      const char * outlevel = ::getenv("PULSAR_OUTPUT_LEVEL");
      m_OutputLevel = atoi(outlevel);
      if ((m_OutputLevel<0) || (m_OutputLevel>2))
	m_OutputLevel=1;
    }
  else
    {
      m_OutputLevel=1;
    }

  if (::getenv("PULSARDATA"))
    {
      const char * psrdata = ::getenv("PULSARDATA");
      m_pulsardata_dir = std::string(psrdata);
    }
  else if (::getenv("SKYMODEL_DIR"))
    {
      const char * psrdata = ::getenv("SKYMODEL_DIR");
      m_pulsardata_dir =  std::string(psrdata)+"/pulsars";
    }
  else
    {
      const char * psrdata = ::getenv("PULSARROOT");
      m_pulsardata_dir =  std::string(psrdata)+"/data";
    }

  if (m_OutputLevel>1)
    {
      WriteToLog("PULSARDATA used is: "+std::string(m_pulsardata_dir));
    }

  // Determine start and end time for FT2 or orbit file 
  try {
    const astro::PointingHistory & history = astro::GPS::instance()->history();
    //    throw(NoHistoryError);
    m_FT2_startMET = history.startTime();
    m_FT2_stopMET = history.endTime();
    m_UseFT2 = 1;
  } catch (astro::GPS::NoHistoryError & )
    {
      m_FT2_startMET = Spectrum::startTime();
      m_FT2_stopMET = astro::GPS::instance()->endTime();
    }

  m_Sim_startMET = Spectrum::startTime();
  m_Sim_stopMET = astro::GPS::instance()->endTime();

  //Init SolarSystem stuffs useful for barycentric decorrections
  astro::JulianDate JDStart(StartMissionDateMJD+JDminusMJD);
  m_earthOrbit = new astro::EarthOrbit(JDStart);   
  astro::SkyDir m_PulsarDir(m_RA,m_dec,astro::SkyDir::EQUATORIAL);
  m_PulsarVectDir = m_PulsarDir.dir();
  m_GalDir = std::make_pair(m_PulsarDir.l(),m_PulsarDir.b());
  m_l = m_GalDir.first;
  m_b = m_GalDir.second;

  //Load Pulsar General data from PulsarDataList.txt
  LoadPulsarData(m_pulsardata_dir,0);
  
  if (!(::getenv("PULSAR_NO_DB")))
    {
      //Save the output txt file..
      int DbFlag = saveDbTxtFile();
      if ((DbFlag !=1) && (m_OutputLevel>1))
	WriteToLog("WARNING! Problem in saving SimPulsars_spin.txt");
    }


  //Binary Pulsar Data
  if (m_BinaryFlag ==1)
    {
      LoadPulsarData(m_pulsardata_dir,1); //loading orbital data from PulsarBinDataList.txt
	
      if (!(::getenv("PULSAR_NO_DB")))
	{
	  int BinDbFlag = saveBinDbTxtFile();
	  if ((BinDbFlag !=1) && (m_OutputLevel>1))
	    WriteToLog("WARNING! Problem in saving SimPulsars_bin.txt");
	}
    }

  if (m_TimingNoiseModel != 0)
    { 
      InitTimingNoise();
    }

  //Instantiate an object of PulsarSim class
  m_Pulsar    = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period);
 
  //Instantiate an object of SpectObj class
  if (m_model == 1)
    {
      m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_ppar0), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
      m_spectrum->SetAreaDetector(EventSource::totalArea());
    }
  else 
    if (m_model == 2)
      {
	m_spectrum = new SpectObj(m_Pulsar->PSRShape(m_PSRShapeName,m_ppar0),1);
	m_spectrum->SetAreaDetector(EventSource::totalArea());
      }
    else
      {
	std::cout << "ERROR!  Model choice not implemented " << std::endl;
	exit(1);
      }


  if (m_OutputLevel > 0)
    {
      
      //Redirect output to a subdirectory
      const char * pulsarOutDir = ::getenv("PULSAROUTFILES");
      
      // override obssim if running in Gleam environment
      if( pulsarOutDir!=0) 
	m_LogFileName = std::string(pulsarOutDir) + "/" + m_PSRname + "Log.txt";
      else
	m_LogFileName = m_PSRname + "Log.txt";

      char temp[200];
      sprintf(temp,"**  Output Level set to: %d",m_OutputLevel);
      WriteToLog(std::string(temp));

      WritePulsarLog();


    }



}

/////////////////////////////////////////////////
PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete m_earthOrbit;
  delete m_PSpectrumRandom;
}

/////////////////////////////////////////////////
double PulsarSpectrum::flux(double time) const
{
  double flux;	  
  flux = m_spectrum->flux(time,m_enphmin);
  return flux;
}


/////////////////////////////////////////////////
/*!
 * \param time time to start the computing of the next photon
 *
 * This method find the time when the next photon will come. It takes into account decorrections due to
 * ephemerides and period derivatives and barycentric decorrections. This method also check for 
 * changes in the ephemerides validity range.
 * For the ephemerides decorrections the parameters used are:
 * <ul>
 * <li> <i>Epoch</i>, espressed in MJD;</li>
 * <li> <i>Initial Phase Phi0</i>;</li>
 * <li> <i>First Period derivative</i>;</li>
 * <li> <i>Second Period derivative</i>;</li>
 * </ul>
 * <br>
 * The barycentric decorrections takes into account conversion TDB->TT, Geometrical and Shapiro delays.
 * This method also calls the interval method in SpectObj class to determine the next photon in a rest-frame without 
 * ephemerides effects.
 */
double PulsarSpectrum::interval(double time)
{  

  double timeTildeDecorr = time + (StartMissionDateMJD)*SecsOneDay; //Arrival time decorrected

  //check if time is before FT2 start
  if (time < m_FT2_startMET)
    {
      timeTildeDecorr+=(m_FT2_startMET-time)+510.;
    }

  //this should be corrected before applying barycentryc decorr + ephem de-corrections

  double timeTilde = 0;

  timeTilde = timeTildeDecorr + getBaryCorr(timeTildeDecorr,0); 

  double timeTildeDemodulated =0.;
  
  //binary demodulation
  if (m_BinaryFlag ==0)
    {
      timeTildeDemodulated = timeTilde;
    }
  else
    {
     
      if (::getenv("PULSAR_OUT_BIN"))
	{
	  timeTildeDemodulated = getIterativeDemodulatedTime(timeTilde,1);      
	}
      else
	{
	  timeTildeDemodulated = getIterativeDemodulatedTime(timeTilde,0);      
	}

      //      std::cout <<"\n****Test inverso:"<< timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay << " con correzioni" << getBinaryDemodulationInverse(timeTildeDemodulated)-(StartMissionDateMJD)*SecsOneDay << std::endl;
      //    std::cout << " indietro -->" << getIterativeDemodulatedTime(getBinaryDemodulationInverse(timeTildeDemodulated),0)-(StartMissionDateMJD)*SecsOneDay << std::endl;
    }

  //Phase assignment

  double intPart=0.; //Integer part
  //  double PhaseNoNoise,PhaseWithNoise=0.;

  //Apply timing noise
  if (m_TimingNoiseModel !=0)
    {
      ApplyTimingNoise(timeTildeDemodulated);
    }

  double initTurns = getTurns(timeTildeDemodulated); //Turns made at this time
  double tStart = modf(initTurns,&intPart)*m_period; // Start time for interval
  if (tStart < 0.)
    tStart+=m_period;

  if (DEBUG)
    {
      std::cout << std::setprecision(30) << "\n" << timeTilde -(StartMissionDateMJD)*SecsOneDay 
		<<  " turns are " << getTurns(timeTilde) 
		<<  " phase is " << tStart/m_period << " phi0 is " << m_phi0 << std::endl;
    }

  CheckEphemeridesValidity(timeTildeDemodulated,initTurns);;

  if (DEBUG)
    {
      if ((int(timeTilde - (StartMissionDateMJD)*SecsOneDay) % 1000) < 1.5)
	std::cout << "\n\n** " << m_PSRname << " Time is: " 
		  << timeTildeDecorr-(StartMissionDateMJD)*SecsOneDay << " s MET in TT | "
		  << timeTilde-(StartMissionDateMJD)*SecsOneDay << "s. MET in SSB in TDB (barycorr.) | "
		  << timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay 
		  << "s. MET in SSB in TDB (barycorr. + demod.)| "
		  << std::endl;
    }

  //Log out the barycentric corrections  
  double bl=0;
  if (::getenv("PULSAR_OUT_BARY"))
    bl = getBaryCorr(timeTilde,1);

  //In case of binary add a modulation...
  if (m_BinaryFlag ==1)
    {
      double ModFactor = 0.5;
      double intOrb;
      //      double Modulation = 1+ModFactor*sin(2*M_PI*modf(initTurns,&intOrb)-DegToRad*m_omega+0.5*M_PI);
      double OrbPhase = getOrbitalPhase(timeTildeDemodulated);

      double Modulation = 1+ModFactor*std::sin(2*M_PI*OrbPhase-DegToRad*m_omega+0.5*M_PI);
      if (DEBUG)
	std::cout << std::setprecision(12)<<timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay<< " OrbPhase " << OrbPhase << std::endl;

      //std::cout << std::setprecision(12)<<TrueAn<<" " <<modf(TrueAn,&intOrb)<<" " << Modulation<<std::endl;
      m_spectrum->SetFluxFactor(Modulation);
    }


  //Ephemerides calculations...
  double inte = m_spectrum->interval(tStart,m_enphmin); //deltaT (in system where Pdot = 0

  


  double inteTurns = inte/m_period; // conversion to # of turns
  double totTurns = initTurns + inteTurns + m_phi0; //Total turns at the nextTimeTilde

  //Applying timingnoise	     
  
  double nextTimeTildeDemodulated = retrieveNextTimeTilde(timeTildeDemodulated, totTurns, (ephemCorrTol/m_period));
  
  double nextTimeTilde = 0.;
  double nextTimeTildeDecorr = 0.;

  //inverse of binary demodulation and of barycentric corrections

  if (m_BinaryFlag == 0)
    {
      nextTimeTilde = nextTimeTildeDemodulated;
      nextTimeTildeDecorr = getDecorrectedTime(nextTimeTilde); //Barycentric decorrections
    }
  else 
    {
      nextTimeTilde = getBinaryDemodulationInverse(nextTimeTildeDemodulated);
      nextTimeTildeDecorr = getDecorrectedTime(nextTimeTilde); //Barycentric decorrections
    }
 
  if (DEBUG)
    { 
      std::cout << std::setprecision(15) << "\nTimeTildeDecorr at Spacecraft (TT) is: " 
		<< timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay << "s." << std::endl;
      std::cout << "Arrival time after start of the simulation : " 
		<< timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay -Spectrum::startTime() << " s." << std::endl; 
      std::cout << std::setprecision(15) <<"  TimeTilde at SSB (in TDB) is: " 
		<< timeTilde - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
      std::cout << std::setprecision(15) <<"  TimeTildeDemodulated: " 
		<< timeTildeDemodulated - (StartMissionDateMJD)*SecsOneDay 
		<< "sec.; phase=" << modf(initTurns,&intPart) << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTildeDemodulated at SSB (in TDB) is: " 
		<< nextTimeTildeDemodulated - (StartMissionDateMJD)*SecsOneDay << " s." << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTilde at SSB (in TDB) is: " 
		<< nextTimeTilde - (StartMissionDateMJD)*SecsOneDay << " s." << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTilde decorrected (TT) is: " 
		<< nextTimeTildeDecorr - (StartMissionDateMJD)*SecsOneDay << " s." << std::endl;
      std::cout << " corrected is " 
		<< nextTimeTildeDecorr + getBaryCorr(nextTimeTildeDecorr,0) - (StartMissionDateMJD)*SecsOneDay << std::endl;
      std::cout << std::setprecision(15) <<"  -->Interval is " <<  nextTimeTildeDecorr - timeTildeDecorr << std::endl;
    }

  //  double interv = nextTimeTildeDecorr - timeTildeDecorr;

  double interv = nextTimeTildeDecorr-(time + (StartMissionDateMJD)*SecsOneDay);

  if (interv <= 0.)
    interv = m_periodVect[0]/100.;
  

  return interv;
  
}

/////////////////////////////////////////////
/*!
 * \param time time where the number of turns is computed
 *
 * This method compute the number of turns made by the pulsar according to the ephemerides.
 * <br>
 * <blockquote>
 * f(t) = f<sub>0</sub> + f<sub>1</sub>(t - t<sub>0</sub>) +
 *      f<sub>2</sub>/2(t - t<sub>0</sub>)<sup>2</sup>.
 * <br>The number of turns is related to the ephemerides as:<br>
 * dN(t) = N<sub>0</sub> + phi<sub>0</sub> + f<sub>0</sub> + f<sub>1</sub>(t-t<sub>0</sub>) + 1/2f<sub>2</sub>(t-t<sub>0</sub>)<sup>2</sup>
 * </blockquote>
 * <br>where 
 * <blockquote>
 *<ul>
 * <li> t<sub>0</sub> is an ephemeris epoch.In this simulator epoch must be expressed in MJD;</li>
 * <li> N<sub>0</sub> is the number of turns at epoch t<sub>0</sub>;
 * <li> phi<sub>0</sub> is the phase shift at epoch t<sub>0</sub>;
 * <li>f<sub>0</sub> pulse frequency at time t<sub>0</sub> (usually given in Hz)</li>,
 * <li>f<sub>1</sub> the 1st time derivative of frequency at time t< sub>0</sub> (Hz s<sup>-1</sup>);</li>
 * <li>f<sub>2</sub> the 2nd time derivative of frequency at time t<sub>0</sub> (Hz s<sup>-2</sup>).</li>
 * </ul>
 * </blockquote>
 */
double PulsarSpectrum::getTurns( double time )
{

  double dt = time - m_t0*SecsOneDay;
  return m_phi0 + m_N0 + m_f0*dt + 0.5*m_f1*dt*dt + ((m_f2*dt*dt*dt)/6.0);
}

/////////////////////////////////////////////
/*!
 * \param tTilde Time from where retrieve the nextTime;
 * \param totalTurns Number of turns completed by the pulsar at nextTime;
 * \param err Phase tolerance between totalTurns and the number of turns at nextTime
 *
 * <br>In this method a recursive way is used to find the <i>nextTime</i>. Starting at <i>tTilde</i>, the method returns 
 * nextTime, the time where the number of turns completed by the pulsar is equal to totalTurns (within the choosen tolerance).  
 */
double PulsarSpectrum::retrieveNextTimeTilde( double tTilde, double totalTurns, double err )
{

  double tTildeUp,NTup = 0.;
  double tTildeDown,NTdown = 0;  

  tTildeUp = tTilde;
  tTildeDown = tTilde;  

  NTdown = totalTurns - getTurns(tTildeDown); 
  NTup = totalTurns - getTurns(tTildeUp);


  int u = 0;
  int nIterations = 0;
  while ((NTdown*NTup)>0)
    {
      u++;
      tTildeUp = tTilde+m_period*pow(2.0,u); 
      NTdown = totalTurns - getTurns(tTildeDown); 
      NTup = totalTurns - getTurns(tTildeUp);
    }

  double tTildeMid = (tTildeDown+tTildeUp)/2.0;
  double NTmid = totalTurns - getTurns(tTildeMid);
  
  while(fabs(NTmid) > err )
    { 
      nIterations++;
      if ((fabs(NTmid) < err) || (nIterations > 200)) break;
      NTmid = totalTurns - getTurns(tTildeMid);
      NTdown = totalTurns - getTurns(tTildeDown); 
      if ((NTmid*NTdown)>0)
	{
  
	  tTildeDown = tTildeMid;
	  tTildeMid = (tTildeDown+tTildeUp)/2.0;
	} 
      else 
	{
	  tTildeUp = tTildeMid;
	  tTildeMid = (tTildeDown+tTildeUp)/2.0;
	}
    }

  if (DEBUG)
    {
      std::cout << "**  RetrieveNextTimeTilde " << std::endl;
      std::cout << std::setprecision(30) << "  Stop up is " << tTildeUp << " NT " << NTup << std::endl;
      std::cout << std::setprecision(30) << "        down is " << tTildeDown << " NT " << NTdown << " u= " << u << std::endl;
      std::cout << std::setprecision(30) << "        mid is " << tTildeMid << " NT " << NTmid << " u= " << u << std::endl;
      std::cout << "     nextTimeTilde is " << tTildeMid << std::endl;
    }

  return tTildeMid;
}

/////////////////////////////////////////////
/*!
 * \param ttInput Photon arrival time at spacecraft in Terrestrial Time (expressed in MJD converted in seconds)
 *
 * <br>
 * This method computes the barycentric corrections for the photon arrival time at spacecraft and returns 
 * the time in Solar System Barycentric Time. The corrections implemented at the moment are:
 * <ul>
 *  <li> Conversion from TT to TDB;
 *  <li> Geometric Time Delay due to light propagation in Solar System;
 *  <li> Relativistic Shapiro delay;
 * </ul>   
 */
double PulsarSpectrum::getBaryCorr( double ttInput, int LogCorrFlag)
{

  if (((ttInput-(StartMissionDateMJD)*SecsOneDay)+510)>m_FT2_stopMET)
    {
      if (m_OutputLevel>1)
	{
	  char temp[200];
	  sprintf(temp,"WARNING!!!Arrival time after FT2 end. No barycentric correction! at t=%.f",ttInput);
	  WriteToLog(std::string(temp));
	}
      return 0.;
    }


  //Start Date;
  astro::JulianDate ttJD(StartMissionDateMJD+JDminusMJD);
  ttJD = ttJD+(ttInput - (StartMissionDateMJD)*SecsOneDay)/SecsOneDay;
  if (DEBUG)
    {
      std::cout << std::setprecision(30) << "\nBarycentric Corrections for time " << ttJD << " (JD)" << std::endl;
    }

  // Conversion TT to TDB, from JDMissionStart (that glbary uses as MJDref)
  double tdb_min_tt = m_earthOrbit->tdb_minus_tt(ttJD);

  double timeMET = ttInput - (StartMissionDateMJD)*SecsOneDay;
  CLHEP::Hep3Vector scPos;


  //Exception error in case of time not in the range of available position (when using a FT2 file)
  try {
    //    std::cout<<astro::GPS::time()<<std::endl;
    //astro::GPS::update(timeMET);
    //    std::cout<<astro::GPS::time()<<std::endl;
    astro::GPS::instance()->time(timeMET);
    scPos = astro::GPS::instance()->position(timeMET);
  } catch (astro::PointingHistory::TimeRangeError & )
    {
      if (DEBUG)
	{
	  std::cout << "WARNING! FT2 Time Error at " << std::setprecision(30) << timeMET << std::endl;
	}
  
      if ((timeMET+510)>=m_FT2_stopMET)
	{
	  return 0.;
	}
      
    }
  //Correction due to geometric time delay of light propagation 
  //GLAST position
  CLHEP::Hep3Vector GeomVect = (scPos/clight);
  double GLASTPosGeom = GeomVect.dot(m_PulsarVectDir);
  double GeomCorr = 0;
  double EarthPosGeom =0.;
  try{
    //Earth position
    GeomVect = - m_solSys.getBarycenter(ttJD);
    EarthPosGeom = GeomVect.dot(m_PulsarVectDir);
    GeomVect = (scPos/clight) - m_solSys.getBarycenter(ttJD);
    GeomCorr = GeomVect.dot(m_PulsarVectDir);
  } catch (astro::SolarSystem::BadDate & )
    {
      std::cout << "WARNING! JD Bad Date " << ttJD <<std::endl;
      GeomCorr = 0.;  
    }


  //Correction due to Shapiro delay.
  CLHEP::Hep3Vector sunV = m_solSys.getSolarVector(ttJD);

  // Angle of source-sun-observer
  double costheta = - sunV.dot(m_PulsarVectDir) / ( sunV.mag() * m_PulsarVectDir.mag() );
  double m = 4.9271e-6; // m = G * Msun / c^3
  double ShapiroCorr = 2.0 * m * log(1+costheta);
  if (DEBUG)
    {
      std::cout << std::setprecision(20) << "** --> TDB-TT = " << tdb_min_tt << std::endl;
      std::cout << std::setprecision(20) << "** --> Geom. delay = " << GeomCorr << std::endl;
      std::cout << std::setprecision(20) << "** --> Shapiro delay = " << ShapiroCorr << std::endl;
      std::cout << std::setprecision(20) << "** ====> Total= " << tdb_min_tt + GeomCorr + ShapiroCorr  << " s." <<std::endl;
    }  


  if ((LogCorrFlag == 1) && (ttInput-StartMissionDateMJD*SecsOneDay > 0.))
    {
      double intPart=0.;
      double PhaseOut = getTurns(ttInput);
      PhaseOut = modf(PhaseOut,&intPart); // phase for linear evolution 
      if (PhaseOut < 0)
	PhaseOut++;



	  std::ofstream BaryCorrLogFile((m_PSRname + "BaryCorr.log").c_str(),std::ios::app);
	  BaryCorrLogFile << std::setprecision(20)
			  << timeMET << "\t" <<  PhaseOut << "\t" << GLASTPosGeom << "\t"<< EarthPosGeom << "\t" 
			  << GeomCorr << "\t" << tdb_min_tt << "\t" << ShapiroCorr << std::endl;
	  BaryCorrLogFile.close();
    }

  return tdb_min_tt + GeomCorr + ShapiroCorr; //seconds
  
}

/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * \param LogFlag Flag for writing output
 *
 * <br>
 * This method compute the binary demodulation in an iterative way 
 * using the getBinaryDemodulation method
 */
double PulsarSpectrum::getIterativeDemodulatedTime(double tInput, int LogFlag)
{

  double timeDemodulated = tInput+getBinaryDemodulation(tInput,0);
            
  double TargetTime = tInput;
  double delay = getBinaryDemodulation(tInput,0);

  int i=0;
  while ((fabs(timeDemodulated-delay-TargetTime) > DemodTol) && ( i < 50))
    {
      i++;
      timeDemodulated = TargetTime+delay;
      delay = getBinaryDemodulation(timeDemodulated,0);
      // std::cout << "Step " << i << " t= " << timeDemodulated-(StartMissionDateMJD)*SecsOneDay<<std::endl;
    }
  
  if (LogFlag == 1)
    {
      delay = getBinaryDemodulation(timeDemodulated,1);
    }
  
  return timeDemodulated;
}

/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * <br>
 * This method returns the true anomaly
 */
double PulsarSpectrum::getTrueAnomaly(double tInput)
{
  //double BinaryRoemerDelay = 0.;
  double dt = tInput-m_t0PeriastrMJD*SecsOneDay;
  double OmegaMean = 2*M_PI/m_Porb;
  double EccAnConst = OmegaMean*(dt - 0.5*(dt*dt*(m_Porb_dot/m_Porb)));
  
  //Calculate Eccenctric Anomaly solving Kepler equation using the atKepler function
  double EccentricAnomaly = 0.; 
  int status = atKepler(EccAnConst, m_ecc, &EccentricAnomaly); //AtKepler
  
  // atKepler not converged
  if (0 != status) {
     throw std::runtime_error("atKepler did not converge.");
  } 
  
  //Calculate True Anomaly
  double TrueAnomaly = 2.0 * std::atan(std::sqrt((1.0+m_ecc)/(1.0-m_ecc))*std::tan(EccentricAnomaly*0.5));
  TrueAnomaly = TrueAnomaly + 2*M_PI*floor((EccentricAnomaly - TrueAnomaly)/ (2*M_PI));
  while ((TrueAnomaly - EccentricAnomaly) > M_PI) TrueAnomaly -= 2*M_PI;
  while ((EccentricAnomaly - TrueAnomaly) > M_PI) TrueAnomaly += 2*M_PI;

  return TrueAnomaly;

}


/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * <br>
 * This method returns the true anomaly
 */
double PulsarSpectrum::getOrbitalPhase(double tInput)
{
  //double BinaryRoemerDelay = 0.;
  double dt = tInput-m_t0PeriastrMJD*SecsOneDay;

  double dp = dt/m_Porb;

  // Compute the complete phase.
  double phase = dp * (1. - dp *m_Porb_dot / 2.0);

  double intOrb;
  double OrbPhase = modf(phase,&intOrb);

  if (OrbPhase<0)
    OrbPhase+=1;

  return OrbPhase;
}


/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * \param LogFlag Flag for writing output
 *
 * <br>
 * This method compute the binary demodulation using the orbital parameters
 * of the pulsar. The corrections computed are:
 * 
 *<ul>
 * <li> Roemer delay
 * <li> Einstein delay
 * <li> Shapiro delay
 *</ul>
 */
double PulsarSpectrum::getBinaryDemodulation( double tInput, int LogDemodFlag)
{
  double BinaryRoemerDelay = 0.;
  double dt = tInput-m_t0PeriastrMJD*SecsOneDay;
  double OmegaMean = 2*M_PI/m_Porb;
  double EccAnConst = OmegaMean*(dt - 0.5*(dt*dt*(m_Porb_dot/m_Porb)));
  
  //Calculate Eccenctric Anomaly solving Kepler equation using the atKepler function
  double EccentricAnomaly = 0.; 
  int status = atKepler(EccAnConst, m_ecc, &EccentricAnomaly); //AtKepler
  
  // atKepler not converged
  if (0 != status) {
     throw std::runtime_error("atKepler did not converge.");
  } 
  
  //Calculate True Anomaly
  double TrueAnomaly = 2.0 * std::atan(std::sqrt((1.0+m_ecc)/(1.0-m_ecc))*std::tan(EccentricAnomaly*0.5));
  TrueAnomaly = TrueAnomaly + 2*M_PI*floor((EccentricAnomaly - TrueAnomaly)/ (2*M_PI));
  while ((TrueAnomaly - EccentricAnomaly) > M_PI) TrueAnomaly -= 2*M_PI;
  while ((EccentricAnomaly - TrueAnomaly) > M_PI) TrueAnomaly += 2*M_PI;
  
  // Calculate longitude of periastron using m_omega_dot
  double Omega = DegToRad*m_omega + DegToRad*m_omega_dot*(TrueAnomaly/OmegaMean);
  
  // compute projected semimajor axis using x_dot
  double asini = m_asini + m_xdot*dt;
  
  //Binary Roemer delay
  BinaryRoemerDelay = ((std::cos(EccentricAnomaly)-m_ecc)*std::sin(Omega)
			     + std::sin(EccentricAnomaly)*std::cos(Omega)*std::sqrt(1-m_ecc*m_ecc));
  //Einstein Delay
  double BinaryEinsteinDelay = m_gamma*std::sin(EccentricAnomaly);
 
  //Shapiro binary delay
  double BinaryShapiroDelay = -2.0*m_shapiro_r*log(1.-m_ecc*std::cos(EccentricAnomaly)-m_shapiro_s*BinaryRoemerDelay);

  BinaryRoemerDelay = asini*BinaryRoemerDelay;
     
  if (DEBUG)
    {
      std::cout << "\n**  Binary modulations t=" << std::setprecision(15) << tInput-(StartMissionDateMJD)*SecsOneDay
		<< " dt=" << tInput - m_t0PeriastrMJD*SecsOneDay << std::endl;
      std::cout << "**  Ecc.Anomaly=" << EccentricAnomaly << " deg. True Anomaly=" << TrueAnomaly << " deg." << std::endl;
      std::cout << "**  Omega=" << Omega << " rad. Major Semiaxis " << asini << " light-sec" << std::endl;
      std::cout << "**  --> Binary Roemer Delay=" << BinaryRoemerDelay << " s." << std::endl;
      std::cout << "**  --> Binary Einstein Delay=" << BinaryEinsteinDelay << " s." << std::endl;
      std::cout << "**  --> Binary Shapiro Delay=" << BinaryShapiroDelay << " s." << std::endl;
  }
  
  if ((LogDemodFlag==1) && (tInput-StartMissionDateMJD*SecsOneDay > 0.))
    {
      std::ofstream BinDemodLogFile((m_PSRname + "BinDemod.log").c_str(),std::ios::app);
      BinDemodLogFile << std::setprecision(15) 
		      << tInput-StartMissionDateMJD*SecsOneDay << "\t" << tInput - m_t0PeriastrMJD*SecsOneDay << "\t"
		      << EccentricAnomaly << "\t" << TrueAnomaly << "\t"
		      << Omega << "\t" << m_ecc << "\t" << asini << "\t"
		      << BinaryRoemerDelay << "\t" 
		      << BinaryEinsteinDelay << "\t" 
		      << BinaryShapiroDelay << "\t" << std::endl;
      BinDemodLogFile.close();
    }

  return -1.*(BinaryRoemerDelay+BinaryEinsteinDelay+BinaryShapiroDelay);
}

/////////////////////////////////////////////
/*!
 * \param CorrectedTime Photon arrival time at SSB (TDB expressed in MJD)
 *
 * <br>
 * This method returns the correspondent decorrected time starting from a photon arrival time
 * at the Solar System barycenter expressed in Barycentric Dynamical Time. This function uses the bisection method
 * for inverting barycentric corrections <br>
 * The corrections implemented at the moment are:
 * <ul>
 *  <li> Conversion from TT to TDB;
 *  <li> Geometric Time Delay due to light propagation in Solar System;
 *  <li> Relativistic Shapiro delay;
 * </ul>   
 */
double PulsarSpectrum::getDecorrectedTime(double CorrectedTime)
{
  //  double CorrectedTime = 54.2342.;
  double deltaT = 510.;

  double tcurr = CorrectedTime-deltaT;
  double fcurr = tcurr + getBaryCorr(tcurr,0); // fx(tcurr);
  //  double fcurr_ct = CorrectedTime;

  if (DEBUG)
    {
      std::cout << "Get decorrected time from time t=" << CorrectedTime << std::endl;
    }

  //double deltaStep = 10.;
  int s = 0;
  int SignDirection = 0;

  //std::cout << "\n\n***\ntstart= " << tcurr << " -->fcurr= " << fcurr 
  //    << " -->fcurr_CT= " << fcurr_ct << " inital delta=" << deltaStep << std::endl;

  if (CorrectedTime <= fcurr) // function increasing
    {
      //      std::cout << "Case 1: de-crescent function" << std::endl;
      SignDirection = 1;
    }
  else
    {
      //std::cout << "Case 2: crescent function" << std::endl;
      SignDirection = 2;
    }

  double deltaStep =  pow(10.,(2.-s));

      while ( fabs(CorrectedTime-fcurr) > baryCorrTol)
	{

	  if (SignDirection == 1)
	    {
	      while (fcurr > CorrectedTime)
		{ 
		  tcurr = tcurr+deltaStep;
		  fcurr = tcurr + getBaryCorr(tcurr,0); //fx(tcurr);
		  //      std::cout << std::setprecision(30) << " tcurr=>" << tcurr << " fcurr " << fcurr 
		  //	<< "df=" << CorrectedTime-fcurr << std::endl; 
		}
	    }
	  else
	    while (fcurr < CorrectedTime)
	      { 
		tcurr = tcurr+deltaStep;
		fcurr = tcurr + getBaryCorr(tcurr,0); //fx(tcurr);
		//std::cout << std::setprecision(30) << " tcurr=>" << tcurr << " fcurr " << fcurr 
		//  << "df=" << CorrectedTime-fcurr << std::endl; 
	      }
	    
	  
	  tcurr = tcurr-deltaStep;
	  fcurr = tcurr + getBaryCorr(tcurr,0); //fx(tcurr);
	  s++;
	  deltaStep =  pow(10.,(2.-s));
	  //std::cout << "\n" << m_PSRname << "  Target superated, decreasing to " << deltaStep 
	  //    << " and starting again1 from " << tcurr << std::endl;
	  //std::cout << std::setprecision(30) << CorrectedTime << 
	  //" <--" << fcurr << " df=" << CorrectedTime-fcurr << std::endl;
	  if ((s > 10) || (deltaStep < 1e-7))
	    {
	      //std::cout << "Skipping ! " << std::endl;
	      break;
	    }
	}

      if (DEBUG)
	{
	  std::cout << std::setprecision(30) <<  s << " -> " << CorrectedTime 
		    << " <--" << fcurr 
		    << " df=" << CorrectedTime-fcurr << std::endl;
	}


      return tcurr;
}

/////////////////////////////////////////////
/*!
 * \param CorrectedTime Photon arrival time Modulatedat SSB (TDB expressed in MJD)
 *
 * <br>
 * This method returns the correspondent inverse-demodulated time of each photons, includingRoemer delay, Einstein delay 
 * and Shapiro delay
 */
double PulsarSpectrum::getBinaryDemodulationInverse( double CorrectedTime)
{

  double deltaMin = m_asini*(1.+m_ecc);//-m_asini*std::sqrt(1.-m_ecc*m_ecc);

  double tcurr = CorrectedTime-deltaMin;
  double fcurr = getIterativeDemodulatedTime(tcurr,0);
  //  double fcurr_ct = CorrectedTime;

  if (DEBUG)
    {
      std::cout << "Get inverse demodulated time for t=" << CorrectedTime << std::endl;
    }

  //double deltaStep = 10.;
  int s = 0;
  int SignDirection = 0;
  double StartStep = 1.0*int(log10(fabs(deltaMin/10.)));
  double deltaStep =  pow(10.,(StartStep-s));
  //std::cout << "START " << StartStep << " dmin " << deltaMin << std::endl;
  if (CorrectedTime <= fcurr) // function increasing
    {
      //      std::cout << "Case 1: de-crescent function" << std::endl;
      SignDirection = 1;
    }
  else
    {
      //std::cout << "Case 2: crescent function" << std::endl;
      SignDirection = 2;
    }

  /*
  std::cout << "\n\n***\n"<<m_PSRname<<"tstart= " << tcurr << " -->fcurr= " << fcurr 
	    << " -->fcurr_CT= " << CorrectedTime << " inital delta=" 
	    << deltaStep << " sign" << SignDirection << std::endl;
  */

  //double tModMid=0.;
  int nMaxStepIterations=0;	
      while ( fabs(CorrectedTime-fcurr) > InverseDemodTol)
	{

          nMaxStepIterations=0;
	  if (SignDirection == 1)
	    {
	      while ((fcurr > CorrectedTime) && (nMaxStepIterations <1000))
		{ 
		  tcurr = tcurr+deltaStep;
		  fcurr = getIterativeDemodulatedTime(tcurr,0);
	          nMaxStepIterations++;	 
                  if (fcurr < CorrectedTime)
                   {
		     break;
                   }
                  //std::cout << std::setprecision(30) << nMaxStepIterations <<" tcurr=>" << tcurr << " fcurr " << fcurr 
		  //	<< "df=" << CorrectedTime-fcurr << std::endl; 
		}
	    }
	  else
	    while ((fcurr < CorrectedTime) && (nMaxStepIterations <1000))
	      { 
		tcurr = tcurr+deltaStep;
		fcurr = getIterativeDemodulatedTime(tcurr,0);
                nMaxStepIterations++;
		//std::cout << std::setprecision(30) << nMaxStepIterations << " tcurr=>" << tcurr << " fcurr " << fcurr 
		  //<< "df=" << CorrectedTime-fcurr << std::endl; 
             if (fcurr > CorrectedTime)
                   {
		     break;
                   }
	      

}
	    
	  
	  tcurr = tcurr-deltaStep;
	  fcurr = getIterativeDemodulatedTime(tcurr,0);
	  s++;
	  deltaStep =  pow(10.,(StartStep-s));
	  //std::cout << "\n" << m_PSRname << "  Target superated, decreasing to " << deltaStep 
	    //  << " and starting again1 from " << tcurr << std::endl;
	  //std::cout << std::setprecision(30) << CorrectedTime << 
	  //" <--" << fcurr << " df=" << CorrectedTime-fcurr << std::endl;
	  if ((s > 8) || (deltaStep < 1e-7))
	    {
	      //std::cout << "Skipping ! " << std::endl;
	      break;
	    }
	}

      //      std::cout << std::setprecision(30) <<  s << " -> " << CorrectedTime 
      //	<< " <--" << fcurr 
      //	<< " df=" << CorrectedTime-fcurr << std::endl;

      return tcurr;
}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method gets from the ASCII <i>PulsarDatalist.txt</i> file stored in <i>/data</i> directory)
 * the names of the files that contain the pulsars parameters. The default one is <i>BasicDataList.txt</i>.
 * This method returns a integer status code (1 is Ok, 0 is failure)
 */
int PulsarSpectrum::getPulsarFromDataList(std::string sourceFileName)
{

  double startTime = 0.;//Spectrum::startTime();
  double endTime = 0.;//astro::GPS::instance()->endTime();

  if (m_UseFT2 == 1)
    {
      startTime =  m_FT2_startMET + (550/86400.);
      endTime = m_FT2_stopMET - (550/86400.);
    } else
      {
	startTime =  m_Sim_startMET + (550/86400.);
	endTime = m_Sim_stopMET - (550/86400.);
      }



  int Status = 0;
  std::ifstream PulsarDataTXT;
  
  if (DEBUG)
    {
      std::cout << "\nOpening Pulsar Datalist File : " << sourceFileName << std::endl;
    }
  
  PulsarDataTXT.open(sourceFileName.c_str(), std::ios::in);
  
  if (! PulsarDataTXT.is_open()) 
    {
      throw "Error! Cannot open file ";
    }
  else {
     char aLine[400];  
    PulsarDataTXT.getline(aLine,400); 
    
    char tempName[15] = "";
    
    double flux, ephem0, ephem1, ephem2, t0Init, t0, t0End, txbary, phi0, period, pdot, p2dot, f0, f1, f2, phas;
    char ephemType[15] = "";
    int tnmodel, binflag;
    
    while (PulsarDataTXT.eof() != 1)
      {

	PulsarDataTXT >> tempName >> flux >> ephemType >> ephem0 >> ephem1 >> ephem2 
		      >> t0Init >> t0 >> t0End  >> txbary >> tnmodel >> binflag;


	//	std::cout << tempName << flux << ephemType << ephem0 << ephem1 << ephem2 
	//  << t0Init << t0 << t0End  << txbary << tnmodel << binflag <<std::endl;

	if (std::string(tempName) == m_PSRname)
	  {
	    
	    Status = 1;
	    m_flux = flux;
	    m_TimingNoiseModel = tnmodel;
	    m_ephemType = ephemType;
	    m_BinaryFlag = binflag;
	    
	    //Check if txbary or t0 are before start of the simulation
	    double startMJD = StartMissionDateMJD+(startTime/86400.)+(550/86400.);
	    // std::cout << "T0 " << t0 << " start " << startMJD << std::endl;
	    if ((t0 < startMJD) || (txbary < startMJD))
	      {
		if (m_OutputLevel>1)
		  {
		    char temp[200];
		    sprintf(temp,"Warning! Epoch t0 out the simulation range (t0-tStart=%.10f) s.: changing to MJD=%d",(t0-startMJD),startMJD);
		    WriteToLog(std::string(temp));
		  }
		t0 = startMJD;
		txbary = startMJD;
	      }
	    
	    //Check if txbary or t0 are after start of the simulation
	    double endMJD = StartMissionDateMJD+(endTime/86400.)-(550/86400.);
	    // std::cout << "T0 " << t0 << " start " << startMJD << std::endl;
	    if ((t0 > endMJD) || (txbary > endMJD))
	      {
		if (m_OutputLevel>1)
		  {

		    char temp[200];
		    sprintf(temp,"Warning! Epoch t0 out the simulation range (t0-tEnd=%.10f) s.: changing to MJD=%d",
			    (t0-endMJD),endMJD);
		    WriteToLog(std::string(temp));
		    sprintf(temp,"**  Tnd at %d corresp. to MJD %d",endTime,endMJD);
		    WriteToLog(std::string(temp));
		  }
		t0 = endMJD;
		txbary = endMJD;
	      }
	    
	    m_t0InitVect.push_back(t0Init);
	    m_t0Vect.push_back(t0);
	    m_t0EndVect.push_back(t0End);
	    m_txbaryVect.push_back(txbary);
	    
	    //Period-type ephemerides
	    if (std::string(ephemType) == "P")
	      {
		
		period = ephem0;
		pdot = ephem1;
		p2dot = ephem2;
		f0 = 1.0/period;
		f1 = -pdot/(period*period);
		f2 = 2*pow((pdot/period),2.0)/period - p2dot/(period*period);
		
	      } 
	    else if (std::string(ephemType) == "F")
	      {
		//Frequency-style ephemrides
		f0 = ephem0;
		f1 = ephem1;
		f2 = ephem2;
		period = 1.0/f0;
	      }
	    
	    m_periodVect.push_back(period);
	    m_pdotVect.push_back(pdot);
	    m_p2dotVect.push_back(p2dot);
	    m_f0Vect.push_back(f0);
	    m_f1Vect.push_back(f1);
	    m_f2Vect.push_back(f2);
	    
	    double dt = (txbary-t0)*SecsOneDay;
	    
	    phi0 = -1.0*(f0*dt
			 + (f1/2.0)*dt*dt
			 + (f2/6.0)*dt*dt*dt); 
	    
	    phi0 = modf(phi0,&phas);
	    
	    if (phi0 < 0. ) 
	      phi0++;
	    
	    m_phi0Vect.push_back(phi0);
	    
	  }
      }
  }

  return Status;
}


/////////////////////////////////////////////
/*!
 * \param pulsar_data_dir: $PULSARDATA directory;
 * \param DataType: Type of data to be loaded (0=General data; 1=Orbital Data);
 *
 * <br>
 * This method load pulsar general data (flux, spin parameters, etc..) of
 * orbital data (kepler parameters, etc..) according to DataType parameter
 * Datafiles containing parameters must be specified in $PULSARDATA/PulsarDataList.txt 
 * for general data or $PULSARDATA/PulsarBinDataList.txt for orbital data
 *
 */
void PulsarSpectrum::LoadPulsarData(std::string pulsar_data_dir, int DataType = 0)
{

  //Scan PulsarDataList.txt for Pulsar general data

  std::string ListFileName;

  if (DataType == 0)
    {
      ListFileName = facilities::commonUtilities::joinPath(pulsar_data_dir, "PulsarDataList.txt");
      if (DEBUG)
	{
	  std::cout << "Reading General Data" << std::endl;
	}
    }
  else if (DataType == 1)
    {
      ListFileName = facilities::commonUtilities::joinPath(pulsar_data_dir, "PulsarBinDataList.txt");
      if (DEBUG)
	{      
	  std::cout << "Reading Orbital Data" << std::endl;
	}
    }

  //new block for reading lines

  //max
 try
   {
     CheckFileExistence(ListFileName);
     
     std::ifstream ListFile(ListFileName.c_str());
     std::string dataline;
     std::vector<std::string> datalines;
     while (std::getline(ListFile, dataline, '\n'))
       {
	 if (dataline != "" && dataline != " "
	     && dataline.find_first_of('#') != 0)
	   {
	     datalines.push_back(dataline);
	   }
       }

     // After lines are parsed, each DataList is processed...
     int PulsarFound=0;
     int l=0;

     while ((l<datalines.size()) && (PulsarFound!=1))
       {
	std::string CompletePathFileName = facilities::commonUtilities::joinPath(pulsar_data_dir,datalines[l]);
	try
	    {

	      if (DataType == 0)
		{
		  PulsarFound = getPulsarFromDataList(CompletePathFileName);
		}
	      else if (DataType == 1)
		{
		  PulsarFound = getOrbitalDataFromBinDataList(CompletePathFileName);
		}

	    }
	 catch (char const *error) //DataList does not exists....
	   {
	     if (m_OutputLevel>1)
	       {
		 WriteToLog("WARNING!"+std::string(error)+std::string(CompletePathFileName)
			    +"... skip to next ASCII Data list file");
	       }
	   }
	
	 l++;
       }

     if (PulsarFound == 0)//if no Datalist contains pulsars...
       {
	 std::cout << "ERROR! Unable to find pulsar " << m_PSRname 
		   << " in any DataList. Check $PULSARDATA" << std::endl;
	 exit(1);
       }


     //if requested, initialize log output for barycentric corrections
     if (::getenv("PULSAR_OUT_BARY"))
       {
	 std::ofstream BaryCorrLogFile((m_PSRname + "BaryCorr.log").c_str());
	 BaryCorrLogFile << "\nPulsar" << m_PSRname 
			 << " : Barycentric corrections log generated by PulsarSpectrum" 
			 << std::endl;
	 BaryCorrLogFile << "tMET\tTDB-TT\tGeomDelay\tShapiroDelay" << std::endl;
	 BaryCorrLogFile.close();
       }


     //if requested, initialize log output for barycentric corrections
      if ((m_BinaryFlag ==1) && (::getenv("PULSAR_OUT_BIN")))
	{
	  std::ofstream BinDemodLogFile((m_PSRname + "BinDemod.log").c_str());
	  BinDemodLogFile << "tMET\tdt\tE\tAt\tOmega\tecc\tasini\tdt_Roemer\tdt_einstein\tdt_shapiro" 
			  << std::endl;
	  BinDemodLogFile.close();
	}
     
   }
 catch (char const *error) //DataList.txt of BinDataList.txt does not exists...
    {
      std::cerr << error << ListFileName << std::endl;
      exit(1);
    }

 if (DataType == 0)
   {
     //Assign as starting ephemeris the first entry of the vectors... 
     m_t0Init = m_t0InitVect[0];
     m_t0 = m_t0Vect[0];
     m_t0End = m_t0EndVect[0];
     m_period = m_periodVect[0];
     m_pdot = m_pdotVect[0];
     m_p2dot = m_p2dotVect[0];
     m_f0 = m_f0Vect[0];
     m_f1 = m_f1Vect[0];
     m_f2 = m_f2Vect[0];
     m_f0NoNoise = m_f0Vect[0];
     m_f1NoNoise = m_f1Vect[0];
     m_f2NoNoise = m_f2Vect[0];
     m_phi0 = m_phi0Vect[0];
   }


}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method gets the orbital parameters of the binary pulsar from a 
 * <i>BinDataList/i> file among those specified in the<i>/data/PulsarBinDataList</i> 
 *
 * The orbital parameters are:
 *
 * m_Porb Orbital period;
 * m_asini Projected major semiaxis;
 * m_eccentr eccentricity;
 * m_omega longitude of periastron;
 * m_t0PeriastronMJD  epoch of periastron passage;
 * m_t0AscNodeMJD Epoch of the ascending node;
 *
 * The key for finding pulsar is the name of the files that contain the pulsars parameters. 
 * The default one is <i>BasicDataList.txt</i>.
 * Extra parameters are used to specify a PPN parameterization for General Relativity;
 * This method returns a integer status code (1 is Ok, 0 is failure)
 */
int PulsarSpectrum::getOrbitalDataFromBinDataList(std::string sourceBinFileName)
{
  int Status = 0;
  std::ifstream PulsarBinDataTXT;
  
  if (DEBUG)
    {
      std::cout << "\nOpening Binary Pulsar BinDatalist File : " << sourceBinFileName << std::endl;
    }
  
  PulsarBinDataTXT.open(sourceBinFileName.c_str(), std::ios::in);
  
  if (! PulsarBinDataTXT.is_open()) 
    {
      throw "Error!Cannot open file";
      //      std::cout << "Error opening BinDatalist file " << sourceBinFileName  
      //	<< " (check whether $PULSARDATA is set)" << std::endl; 
      //Status = 0;
      //exit(1);
    }
  
  char aLine[400];  
  PulsarBinDataTXT.getline(aLine,400); 
 
  char tempName[30];
  double porb,asini,ecc,omega,t0peri,t0asc,ppn;

  while (PulsarBinDataTXT.eof() != 1)
    {
      
      PulsarBinDataTXT >> tempName >> porb >> asini >> ecc >> omega >> t0peri >> t0asc >> ppn;
      
      if (std::string(tempName) == m_PSRname)
	{
	  Status = 1;
	  m_Porb = porb;
	  m_asini = asini;
	  m_ecc = ecc;
	  m_omega = omega;
	  m_t0PeriastrMJD = t0peri;
	  m_t0AscNodeMJD = t0asc;
	  m_PPN = ppn;
	}
    }

  return Status;
}


/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * Initialize variable related to timing noise
 */
void PulsarSpectrum::InitTimingNoise()
{

  if (::getenv("PULSAR_OUT_TNOISE"))
    {
      std::ofstream TimingNoiseLogFile((m_PSRname + "TimingNoise.log").c_str());
      TimingNoiseLogFile << "Timing Noise log file using model: " << m_TimingNoiseModel << std::endl;
      TimingNoiseLogFile << "tMET\tA\tS0\tS1\tS2\tf0_l\tf0_n\tf1_l\tf1_n\tf2_l\tf2_n\tPhi_l\tPhi_n"
			 <<std::endl;
      TimingNoiseLogFile.close();
    }
  
  //define a default mean rate of about 1 day
  m_TimingNoiseMeanRate = 1/86400.;
  // according to Poisson statistics
  double startTime = Spectrum::startTime();
  //Determine next Timing noise event according to the rate R m_TimingNoiseRate      
  m_TimingNoiseTimeNextEvent = startTime -log(1.-m_PSpectrumRandom->Uniform(1.0))/m_TimingNoiseMeanRate; 
  
}


/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * TnoiseInputTime input time where to apply noise
 *
 * Do the timing noise calculation on ephemerides
 */
void PulsarSpectrum::ApplyTimingNoise(double TnoiseInputTime)
{

  double intPart=0.; //Integer part
  double PhaseNoNoise,PhaseWithNoise=0.;

  //Check for a next Timing noise event
      if ((TnoiseInputTime-(StartMissionDateMJD)*SecsOneDay) > m_TimingNoiseTimeNextEvent)
	{
	  //interval according to Poisson statistics
	  m_TimingNoiseTimeNextEvent+= -log(1.-m_PSpectrumRandom->Uniform(1.0))/m_TimingNoiseMeanRate;
 
	  if (DEBUG)
	    {
	      std::cout << std::setprecision(30)<<"Timing Noise Event!Next Event at t=" 
			<< m_TimingNoiseTimeNextEvent << " |dt=" 
			<< m_TimingNoiseTimeNextEvent
		-(TnoiseInputTime-(StartMissionDateMJD)*SecsOneDay) <<std::endl;
	    }


	  //Timing noise managing...

	  double S0=0.;
	  double S1=0.;
	  double S2=0.;

	  if (m_TimingNoiseModel ==1) // Timing Noise #1
	    {
	      m_f2 = 0;
	      PhaseNoNoise = getTurns(TnoiseInputTime);
	      PhaseNoNoise = modf(PhaseNoNoise,&intPart); // phase for linear evolution 
	      if ( PhaseNoNoise <0.)
		PhaseNoNoise+=1.;
	      m_TimingNoiseActivity = 6.6 + 0.6*log10(m_pdot) + m_PSpectrumRandom->Gaus(0,0.5);
	      
	      //estimate an f2
	      double Sign = m_PSpectrumRandom->Uniform();
	      if (Sign > 0.5)
		m_f2 = ((m_f0*6.*std::pow(10.,m_TimingNoiseActivity))*1e-24);
	      else
		m_f2 = -1.0*((m_f0*6.*std::pow(10.,m_TimingNoiseActivity))*1e-24);
	    }
	  else if ((m_TimingNoiseModel >1) && (m_TimingNoiseModel < 5)) // Timing Noise RW -Cordes-Downs
	    {

	      m_f2 = 0.;
	      double tempPhi0 = m_phi0; 
	      m_phi0 = 0.;
	      double tempF0 = m_f0;
	      m_f0 = m_f0NoNoise;	      
	      double tempF1 = m_f1;
	      m_f1 = m_f1NoNoise;	      


	      PhaseNoNoise = getTurns(TnoiseInputTime);
	      PhaseNoNoise = modf(PhaseNoNoise,&intPart); // phase for linear evolution 
	      if ( PhaseNoNoise <0.)
		PhaseNoNoise+=1.;

	      //Determine Crab RMS
	      double dt_days = (TnoiseInputTime-(StartMissionDateMJD)*SecsOneDay)/SecsOneDay;
	      double s_rms_crab = 0.012*pow((dt_days/1628),1.5);

	      //Activity parameter
	      //m_TimingNoiseActivity = -1.37+0.71*log10(m_pdot*1E15);
	      //double s_rms = pow(10.0,m_TimingNoiseActivity)*s_rms_crab;
              m_TimingNoiseActivity = -1.37+0.71*log(m_pdot*1E15);
              double s_rms = exp(m_TimingNoiseActivity)*s_rms_crab;

	      if (m_TimingNoiseModel ==2) //Case 1 :PN
		{
		  S0 = (3.7*3.7*s_rms*s_rms)*(2/(SecsOneDay*dt_days));
		  m_TimingNoiseRMS = std::sqrt((S0/m_TimingNoiseMeanRate));
		  m_phi0 = tempPhi0+m_PSpectrumRandom->Gaus(0,m_TimingNoiseRMS);
		  // std::cout << std::setprecision(30) << "A="<<m_TimingNoiseActivity<<" dt"<<dt_days << " s_rms" << s_rms << 
		  //" S0 " << S0 << " RMS " << m_TimingNoiseRMS << " phi0 " << m_phi0 << std::endl;
		}
	      else if (m_TimingNoiseModel ==3) //Case 2 :PN
		{
		  S1 =(15.5*15.5*s_rms*s_rms)*(12./pow((SecsOneDay*dt_days),3));
		  m_TimingNoiseRMS = std::sqrt((S1/m_TimingNoiseMeanRate));
		  m_f0 = tempF0+m_PSpectrumRandom->Gaus(0,m_TimingNoiseRMS);
		  //std::cout << std::setprecision(30) << "A="<<m_TimingNoiseActivity<<" dt"<<dt_days << " s_rms" << s_rms << 
		  //" S1 " << S1 << " RMS " << m_TimingNoiseRMS << " f0 " << m_f0 << std::endl;
		}
	      else if (m_TimingNoiseModel ==4) //Case 3 :SN
		{
		  S2 =(23.7*23.7*s_rms*s_rms)*(120./pow((SecsOneDay*dt_days),5));
		  m_TimingNoiseRMS = std::sqrt((S2/m_TimingNoiseMeanRate));
		  m_f1 = tempF1+m_PSpectrumRandom->Gaus(0,m_TimingNoiseRMS);
		  //std::cout << std::setprecision(30) << "A="<<m_TimingNoiseActivity<<" dt"<<dt_days << " s_rms" << s_rms << 
		  //  " S2 " << S2 << " RMS " << m_TimingNoiseRMS << " f1 " << m_f1 << std::endl;
		}



	    }

	  if ((DEBUG) || (::getenv("PULSAR_OUT_TNOISE")))
	    {
	      PhaseWithNoise = getTurns(TnoiseInputTime);
	      PhaseWithNoise = modf(PhaseWithNoise,&intPart); // phase for linear evolution 
	      if ( PhaseWithNoise <0.)
		PhaseWithNoise+=1.;
	    }

	  if (::getenv("PULSAR_OUT_TNOISE"))
	    {

	      std::ofstream TimingNoiseLogFile((m_PSRname + "TimingNoise.log").c_str(),std::ios::app);
	      m_f2NoNoise = 0.;
	      double ft_l = GetFt(TnoiseInputTime,m_f0NoNoise,m_f1NoNoise,m_f2NoNoise);
	      double ft_n = GetFt(TnoiseInputTime,m_f0,m_f1,m_f2);
	      double ft1_l = GetF1t(TnoiseInputTime,m_f1NoNoise,m_f2NoNoise);
	      double ft1_n = GetF1t(TnoiseInputTime,m_f1,m_f2);
	      double ft2_l = m_f2NoNoise;//
	      double ft2_n = m_f2;//
	      
	      TimingNoiseLogFile << std::setprecision(30) << TnoiseInputTime-(StartMissionDateMJD)*SecsOneDay
				 << "\t" << m_TimingNoiseActivity
				 << "\t" << S0 
				 << "\t" << S1
				 << "\t" << S2 
				 << "\t" <<ft_l << "\t" << ft_n 
				 << "\t" <<ft1_l << "\t" << ft1_n 
				 << "\t" <<ft2_l << "\t" << ft2_n 			 
				 << "\t" << PhaseNoNoise << "\t" << PhaseWithNoise << std::endl;
	    }
	  

	  if (DEBUG)
	    {

	      std::cout << std::setprecision(30) << " Activity=" << m_TimingNoiseActivity 
			<< "f2=" << m_f2 << " PN=" << PhaseWithNoise 
			<< " dPhi=" << PhaseWithNoise-PhaseNoNoise <<std::endl;
	    }
	}
}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method check if the current time is within valid ephemerides
 */
void PulsarSpectrum::CheckEphemeridesValidity(double EphCheckTime, double initTurns)
{
  //Checks whether ephemerides (period,pdot, etc..) are within validity ranges
  if (((EphCheckTime/SecsOneDay) < m_t0Init) || ((EphCheckTime/SecsOneDay) > m_t0End)) 
    {

      if (m_OutputLevel>1)
	{
	  WriteToLog("WARNING! Time is out of range of validity for pulsar "+std::string(m_PSRname) 
		     +": Switching to new ephemerides set...");
	}


	for (unsigned int e=0; e < m_t0Vect.size();e++)
	  if (((EphCheckTime/SecsOneDay) > m_t0InitVect[e]) && ((EphCheckTime/SecsOneDay) < m_t0EndVect[e])) 
	    {
	
	      m_t0Init = m_t0InitVect[e];
	      m_t0 = m_t0Vect[e];
 	      m_t0End = m_t0EndVect[e];
	      m_f0 = m_f0Vect[e];
	      m_f1 = m_f1Vect[e];
	      m_f2 = m_f2Vect[e];
	      m_f0NoNoise = m_f0Vect[e];
	      m_f1NoNoise = m_f1Vect[e];
	      m_f2NoNoise = m_f2Vect[e];
	      m_period = m_periodVect[e];
	      m_pdot = m_pdotVect[e];
	      m_p2dot = m_p2dotVect[e];
	      m_phi0 = m_phi0Vect[e];

	      if (m_OutputLevel>1)
		{
		  WriteToLog("Valid Ephemerides set found:");
		  char temp[200];
		  sprintf(temp,"MJD(%d-%d) --> Epoch t0 = MJD %d",m_t0Init,m_t0End,m_t0);
		  WriteToLog(std::string(temp));
		  sprintf(temp,"f0: %.f Hz | f1: %.e Hz/s | f2 %.e Hz/s2 ",m_f0,m_f1,m_f2);
		  WriteToLog(std::string(temp));
		  sprintf(temp,"P0: %.f s | P1: %.e s/s | P2 %.e s/s2 ",m_period,m_pdot,m_p2dot);
		  WriteToLog(std::string(temp));
		}


	      //Re-instantiate PulsarSim and SpectObj
	      delete m_Pulsar;

	      m_Pulsar = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period);

	      if (m_model == 1)
		{

		  delete m_spectrum;
		  m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_ppar0), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
		  m_spectrum->SetAreaDetector(EventSource::totalArea());
	      
		}

	      m_N0 = m_N0 + initTurns - getTurns(EphCheckTime); //Number of turns at next t0
	      double intN0;
	      double N0frac = modf(m_N0,&intN0); // Start time for interval
	      m_N0 = m_N0 - N0frac;
	      if (m_OutputLevel > 1)
		{
		  std::cout << std::setprecision(20) << " Turns now are " << initTurns  
			    << " ; at t0 " << m_N0 << std::endl;	           
		}
	      //m_phi0 = m_phi0 - N0frac;
	       
	      if (DEBUG)
		{
		  std::cout << std::setprecision(20) << " At Next t0 Number of Turns will be: " << m_N0 << std::endl;
		}
	    }
	 else
	   {
	     if (m_OutputLevel>1)
	       {
		 WriteToLog("WARNING! Valid ephemerides not found!Proceeding with the old ones");
	       }
	   }
    }
}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method saves the relevant information in a file, named SimPulsar_spin.txt.
 * The format of this ASCII file is such that it can be given to gtpulsardb to produce
 * a D4-compatible FITS file, that can be used with pulsePhase
 */
int PulsarSpectrum::saveDbTxtFile()
{
  int Flag = 0;

  std::string DbOutputFileName = "SimPulsars_spin.txt";
  std::ifstream DbInputFile;

  //Checks if the file exists

  DbInputFile.open(DbOutputFileName.c_str(), std::ios::in);

  if (! DbInputFile.is_open()) 
    {
      Flag = 0;
    }
  else
    {
      Flag = 1;
    }

  DbInputFile.close();

  if (DEBUG)
    {
      std::cout << "Saving Pulsar ephemerides on file " << DbOutputFileName << std::endl;
    }

  std::ofstream DbOutputFile;

  if (Flag == 0)
    {
      DbOutputFile.open(DbOutputFileName.c_str(), std::ios::out);
      DbOutputFile << "# Simulated pulsars output file generated by PulsarSpectrum." << std::endl;
      DbOutputFile << "SPIN_PARAMETERS\n";
      DbOutputFile << "EPHSTYLE = FREQ\n\n# Then, a column header."  << std::endl;
      DbOutputFile << "PSRNAME RA DEC EPOCH_INT EPOCH_FRAC TOAGEO_INT TOAGEO_FRAC TOABARY_INT TOABARY_FRAC ";
      DbOutputFile << "F0 F1 F2 RMS VALID_SINCE VALID_UNTIL BINARY_FLAG SOLAR_SYSTEM_EPHEMERIS OBSERVER_CODE" <<std::endl;
      DbOutputFile.close();
  }

  //Writes out the infos of the file
  DbOutputFile.open(DbOutputFileName.c_str(),std::ios::app);
  double tempInt, tempFract;
  for (unsigned int ep = 0; ep < m_periodVect.size(); ep++)
    {
      DbOutputFile << "\"" << m_PSRname << std::setprecision(10) << "\" " << m_RA << " " << m_dec << " ";
      tempFract = modf(m_t0Vect[ep],&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";
    
      tempFract = getDecorrectedTime(m_txbaryVect[ep]*SecsOneDay)/SecsOneDay;


      tempFract = modf(tempFract,&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";

      tempFract = modf(m_txbaryVect[ep],&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";

      std::string BinFlag;
      if (m_BinaryFlag==0)
	BinFlag = "F";
      else if (m_BinaryFlag==1)
	BinFlag = "T";


      std::string SolarEph;
    
      if (::getenv("PULSAR_EPH"))
	{
	  const char * soleph = ::getenv("PULSAR_EPH");
	  SolarEph = " \"" + std::string(soleph)  + "\"";
	}
      else
	{
	  SolarEph = " \"JPL DE405\"";
	}

      DbOutputFile << std::setprecision(14) 
		   << m_f0Vect[ep] << " " << m_f1Vect[ep] << " " << m_f2Vect[ep] << " " 
		   << m_TimingNoiseRMS*1e3 << " "  
		   << m_t0InitVect[ep] << " " << m_t0EndVect[ep] << " "
		   << BinFlag << SolarEph << " MR" << std::endl;

    }
  

  DbOutputFile.close();

  if (DEBUG) 
    if (Flag == 0)
      { 
	std::cout << "Database Output file created from scratch " << std::endl;
      } 
    else 
      {
	std::cout << "Appendended data to existing Database output file" << std::endl;
      }

  return Flag;
}

/////////////////////////////////////////////
/*!
 * \param NameFileToCheck
 *
 * <br>
 * This method simply check the file and return an exception
 */
void PulsarSpectrum::CheckFileExistence(std::string NameFileToCheck)
{

  std::ifstream FileToCheck;

  FileToCheck.open(NameFileToCheck.c_str(), std::ios::in);

  if (!FileToCheck.is_open()) 
    {
      throw "Error!Cannot open file";
    }

  FileToCheck.close();

}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method saves the relevant orbital parameters in a file, named SimPulsar_bin.txt.
 * The format of this ASCII file is such that it can be given to gtpulsardb to produce
 * a D4-compatible FITS file, that can be used with pulsePhase
 */
int PulsarSpectrum::saveBinDbTxtFile()
{
  int Flag = 0;

  std::string DbBinOutputFileName = "SimPulsars_bin.txt";
  std::ifstream DbBinInputFile;

  //Checks if the file exists

  DbBinInputFile.open(DbBinOutputFileName.c_str(), std::ios::in);

  if (! DbBinInputFile.is_open()) 
    {
      Flag = 0;
    }
  else
    {
      Flag = 1;
    }

  DbBinInputFile.close();

  if (DEBUG)
    {
      std::cout << "Saving Pulsar Orbital Data on file " << DbBinOutputFileName << std::endl;
    }

  std::ofstream DbBinOutputFile;

  if (Flag == 0)
    {
      DbBinOutputFile.open(DbBinOutputFileName.c_str(), std::ios::out);
      DbBinOutputFile << "# Simulated pulsars orbital data output file generated by PulsarSpectrum." << std::endl;
      DbBinOutputFile << "ORBITAL_PARAMETERS\nEPHSTYLE = DD /Simplified model\n# This file can be converted to a D4 fits file using:"  << std::endl;
      DbBinOutputFile << "# >gtpulsardb SimPulsars_bin.txt" << std::endl;
      DbBinOutputFile << "PSRNAME PB PBDOT A1 XDOT ECC ECCDOT OM OMDOT T0 GAMMA SHAPIRO_R SHAPIRO_S OBSERVER_CODE SOLAR_SYSTEM_EPHEMERIS" << std::endl;;
      DbBinOutputFile.close();
  }

  //Writes out the infos of the file
  DbBinOutputFile.open(DbBinOutputFileName.c_str(),std::ios::app);
  
  DbBinOutputFile << "\"" << m_PSRname << "\" ";                                   // pulsar name
  DbBinOutputFile << std::setprecision(10) << m_Porb << " " << m_Porb_dot << " ";   // Orbital period and derivative
  DbBinOutputFile << std::setprecision(10) << m_asini << " " << m_xdot << " ";      // Projected semi-mayor axis and derivative
  DbBinOutputFile << std::setprecision(10) << m_ecc << " " << m_ecc_dot << " ";     // Eccentricity and derivative
  DbBinOutputFile << std::setprecision(10) << m_omega << " " << m_omega_dot << " "; // Long. Periastron and derivative
  DbBinOutputFile << std::setprecision(10) << m_t0PeriastrMJD << " " <<  m_gamma << " "; // t0 of Periastron and PPN gamma
//Shapiro Parameters   
  DbBinOutputFile << std::setprecision(10) << m_shapiro_r << " " << m_shapiro_s; //Shapiro Parameters
  std::string SolarEph;
  if (::getenv("PULSAR_EPH"))
    {
      const char * soleph = ::getenv("PULSAR_EPH");
      SolarEph = " \"" + std::string(soleph)  + "\"";
    }
  else
    {
      SolarEph = " \"JPL DE405\"";
    }

  DbBinOutputFile << " MR "<<SolarEph<<std::endl;// Observer code and ephemerides
 

  DbBinOutputFile.close();


  //In this case a summary D4 file is created
  std::ofstream DbSumInputFile("SimPulsars_summary.txt");
  DbSumInputFile << "SimPulsars_spin.txt\nSimPulsars_bin.txt" <<std::endl;
  DbSumInputFile.close();

  if (DEBUG) 
    if (Flag == 0)
      { 
	std::cout << "Database for Binary pulsars file created from scratch " << std::endl;
      } 
    else 
      {
	std::cout << "Database for Binary pulsars appended to existing binary Database output file" << std::endl;
      }

  return Flag;

}


/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method saves the relevant information about pulsar in a log file
 */
void PulsarSpectrum::WritePulsarLog()
{

  //Write infos to Log file  
  
  std::ofstream PulsarLog(m_LogFileName.c_str(),std::ios::app);

  PulsarLog << "** PulsarSpectrum: " << "********   PulsarSpectrum Log for pulsar" << m_PSRname << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**   Name : " << m_PSRname << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**   Position : (RA,Dec)=(" << m_RA << "," << m_dec 
	    << ") ; (l,b)=(" << m_l << "," << m_b << ")" << std::endl; 
  PulsarLog << "** PulsarSpectrum: "<< "**   Flux above 100 MeV : " << m_flux << " ph/cm2/s " << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**   Enphmin: " << m_enphmin << " keV | Enphmax: " << m_enphmax << " keV" << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;
  

  //Write info about FT2 
  if (m_UseFT2 == 0)
    {
      PulsarLog << "** PulsarSpectrum: "<< "** No FT2 file used " << std::endl;
    }
  else
    PulsarLog << "** PulsarSpectrum: "<< "** FT2 file used " << std::endl;
  
  PulsarLog << "** PulsarSpectrum: "<< "** Start time:" << std::setprecision(30) << m_FT2_startMET << " s. MET |  End time:" 
	    << m_FT2_stopMET << " s. MET" << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "** Simulation start at " << m_Sim_startMET << " s. MET and ends at :" 
	    << m_Sim_stopMET << std::endl;
  PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;

  //Writes down on Log all the ephemerides
  for (unsigned int n=0; n < m_t0Vect.size(); n++)
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Ephemerides valid from " << m_t0InitVect[n] 
		<< " to " << m_t0EndVect[n] << " (MJD): " << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Epoch (MJD) :  " << m_t0Vect[n] << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< std::setprecision(8) << "**     TxBary (MJD) where fiducial point (phi=0) is reached : " 
		<< m_txbaryVect[n] << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Phi0 (at Epoch t0) : " << m_phi0Vect[n] << std::endl;
      
      if (m_ephemType == "P")
	{
	  PulsarLog << "** PulsarSpectrum: "<< "**     Ephemerides type: PERIOD" << std::endl;
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     Period : " 
		    << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     Pdot : " 
		    <<  m_pdotVect[n]  << " | f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     P2dot : " 
		    <<  m_p2dotVect[n]  << " | f2: " << m_f2Vect[n] << std::endl; 
	} 
      else if (m_ephemType == "F")
	{
	  PulsarLog << "** PulsarSpectrum: "<< "**Ephemerides type: FREQUENCY" << std::endl;
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     Period : " 
		    << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << "** PulsarSpectrum: "<< std::setprecision(14) << "**     f2: " << m_f2Vect[n] << std::endl; 
	}
    }

  //MJDRef
  PulsarLog << "** PulsarSpectrum: "<< "**   Mission Reference time: MJD " << StartMissionDateMJD << " (" 
	    << std::setprecision(12) << (StartMissionDateMJD+JDminusMJD)*SecsOneDay 
	    << " sec.)" << std::endl;

  //SimulationModel
  if (m_model ==1)
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Model chosen : " << m_model 
		<< " --> Using Phenomenological Pulsar Model " << std::endl;  
    } else if (m_model == 2)
      {
	PulsarLog << "** PulsarSpectrum: "<< "**   Model chosen : " << m_model 
		  << " --> Using External 2-D Pulsar Shape" << std::endl;  
      }

  PulsarLog << "** PulsarSpectrum: "<< "**   Effective Area set to : " << m_spectrum->GetAreaDetector() << " m^2 " << std::endl; 
  PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;

  //Timing noise
  if (m_TimingNoiseModel == 1) // Timing model #1 - Delta8 parameter (Arzoumanian94)
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Timing Noise Model : " << m_TimingNoiseModel 
		<< " (Stability parameter, Arzoumanian 1994)" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
    }
  else if (m_TimingNoiseModel == 2) //Timing model #2 - PN Random Walk (Cordes-Downs 1985) 
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Timing Noise Model : " << m_TimingNoiseModel 
		<< " (PN Random Walk; Cordes-Downs 1985)" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
    }
  else if (m_TimingNoiseModel == 3) //Timing model #3 - FN Random Walk (Cordes-Downs 1985) 
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Timing Noise Model : " << m_TimingNoiseModel 
		<< " (FN Random Walk; Cordes-Downs 1985)" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
    }
  else if (m_TimingNoiseModel == 4) //Timing model #4 - SN Random Walk (Cordes-Downs 1985) 
    {
      PulsarLog << "** PulsarSpectrum: "<< "**   Timing Noise Model : " << m_TimingNoiseModel 
		<< " (SN Random Walk; Cordes-Downs 1985)" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
    }
  

  //Orbital info
  if (m_BinaryFlag == 1)
    {
      PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**   Pulsar in a Binary System! Orbital Data:" << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Orbital period: " << m_Porb << " s." << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Projected major semiaxis (a * sini): " << m_asini << " lightsec." <<std::endl; 
      PulsarLog << "** PulsarSpectrum: "<< "**     Eccentricity: " << m_ecc << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Longitude of periastron: " <<  m_omega << " deg." << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Epoch of Periastron (MJD): " << m_t0PeriastrMJD << std::endl;
      PulsarLog << "** PulsarSpectrum: "<< "**     Epoch of Ascending Node (MJD): " << m_t0AscNodeMJD << std::endl;
      if (m_PPN ==0)
	PulsarLog << "** PulsarSpectrum: "<< "**   No Post Newtonian Parameterization " << std::endl; 

      PulsarLog << "** PulsarSpectrum: "<< "**************************************************" << std::endl;
    }

  PulsarLog.close();
}

//////////////////////////////////////////////////////////
// no longer used
/////////////////////////////////////////////////////////
double PulsarSpectrum::GetEccentricAnomaly(double mytime)
{
  double OmegaMean = 2*M_PI/m_Porb;
  double dtime = (mytime-m_t0PeriastrMJD*SecsOneDay);
  double EccAnConst = OmegaMean*(dtime - 0.5*(m_Porb_dot/m_Porb)*dtime*dtime);

  double Edown = 0.;
  double EccAnDown = Edown-(m_ecc*std::sin(Edown))-EccAnConst;

  double Eup = 2*M_PI;
  double EccAnUp = Eup-(m_ecc*std::sin(Eup))-EccAnConst;

  double Emid = 0.5*(Eup + Edown);
  double EccAnMid = Emid-(m_ecc*std::sin(Emid))-EccAnConst;

  int i=0;
  while (fabs(EccAnMid) > 5e-7)
    {
      if ((EccAnDown*EccAnMid) < 0)
	{
	  Eup = Emid;
	  EccAnUp = EccAnMid;
	}
      else
	{
	  Edown = Emid;
	  EccAnDown = EccAnMid;
	}
      i++;
      Emid = 0.5*(Eup + Edown);
      EccAnMid = Emid-(m_ecc*std::sin(Emid))-EccAnConst;
      if (fabs(EccAnMid) < 5e-7)
	break;
    }

  return Emid;
}


/////////////////////////////////////////////
double PulsarSpectrum::energy(double time)
{
  return m_spectrum->energy(time,m_enphmin)*1.0e-3; //MeV
}

/////////////////////////////////////////////
/*!
 * \param time input time for calculating f(t)
 *
 * This method computes the frequency at a given instant t
 *
 */ 
double PulsarSpectrum::GetFt(double time, double myf0, double myf1, double myf2)
{
  double dt = time - m_t0*SecsOneDay;
  return myf0 + myf1*dt + 0.5*myf2*dt*dt;
}

/////////////////////////////////////////////
/*!
 * \param time input time for calculating f'(t)
 *
 * This method computes the frequency first derivativeat a given instant t
 *
 */ 
double PulsarSpectrum::GetF1t(double time, double myf1, double myf2)
{
  double dt = time - m_t0*SecsOneDay;
  return myf1 + myf2*dt;
}



/////////////////////////////////////////////
/*!
 * \param input String to be parsed;
 * \param index Position of the parameter to be extracted from input;
 *
 * <br>
 * From a string contained in the XML file a parameter is extracted according to the position 
 * specified in <i>index</i> 
 */
std::string PulsarSpectrum::parseParamList(std::string input, unsigned int index)
{
  std::vector<std::string> output;
  unsigned int i=0;
  int StrLength=input.length();

  while (i<StrLength)
    {
      i=input.find_first_of(",");
      std::string f = ( input.substr(0,i).c_str() );      
      //std::cout << "i=" <<"sub " << f <<std::endl;     
      input=input.substr(i+1);
      output.push_back(f);
    }

  if(index>=output.size()) return "";
  return output[index];

}




/* old code:
std::string PulsarSpectrum::parseParamList(std::string input, unsigned int index)
{
  std::vector<std::string> output;
  unsigned int i=0;
  
  for(;!input.empty() && i!=std::string::npos;){
   
    i=input.find_first_of(",");
    std::string f = ( input.substr(0,i).c_str() );
    output.push_back(f);
    input= input.substr(i+1);
  } 

  if(index>=output.size()) return "";
  return output[index];
}
*/


//////////////////////////////////////////////////
/*!
 * \param Line line to be written in the output log file
 *
 * This method writes a line into a log file
*/
void PulsarSpectrum::WriteToLog(std::string Line)
{
  std::ofstream PulsarLog(m_LogFileName.c_str(),std::ios::app);
  PulsarLog << "** PulsarSpectrum: " << Line << std::endl;
  PulsarLog.close();
}

