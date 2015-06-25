/**
 * @file IRB_routines.cxx
 * @brief
 *
 * Models implemented by Luis C. Reyes... structure inherited from Liz
 * Hays' original code Older EBL models can be found in previous
 * versions of this file.
 *
 * Tau is calculated as a function of energy (GeV) and redshift for
 * the following models in the GLAST range:
 *
 *  EBL model 0: Kneiske et al "Best Fit" (A&A 413, 807-815, (2004))
 *  EBL model 1: Primack, Bullock, Somerville (2005)  astro-ph 0502177
 *  EBL model 2: Kneiske et al "High-UV" (A&A 413, 807, 815, (2004))
 *  EBL model 3: Stecker, Malkan, and Scully (ApJ, 648, 774 (2006) and 
 *               erratum: ApJ 658, 1392S (2007))
 *  EBL model 4: Franceschini et al (2008)  arXiv:0805.1841v2
 *  EBL model 5: Finke et al. (2009) arXiv:0905.1115
 *  EBL model 6: Gilmore et al. (2009)  ArXiv paper, 0905.1144
 *  EBL model 7: Stecker, Malkan, and Scully (ApJ, 648, 774 (2006)  "FAST EVOLUTION"
 *  EBL model 8: Salamon & Stecker (ApJ 1998, 493:547-554) with metallicity correction
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/eblAtten/src/IRB_routines.cxx,v 1.13 2013/12/06 15:20:44 jchiang Exp $
 */


//-------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iostream>

#include "facilities/commonUtilities.h"

#include "AsciiTableModel.h"
#include "Primack05.h"

using namespace std;

namespace IRB {

float calcKneiske(float energy, float redshift);
float calcPrimack05(float energy, float redshift);
float calcKneiske_HighUV(float energy, float redshift);
float calcStecker05(float energy, float redshift);
float calcStecker05_FE(float energy, float redshift);
float calcKneiske_extendedEnergyRange(float energy, float redshift);
float calcKneiske_HighUV_extendedEnergyRange(float energy, float redshift);
float calcFranceschini(float energy, float redshift);
float calcFinke(float energy, float redshift);
float calcGilmore(float energy, float redshift);
float calcSalamonStecker(float energy, float redshift);
float calcGeneric (float energy, float redshift);
float calcGilmore12_fixed(float energy, float redshift);
float calcGilmore12_fiducial(float energy, float redshift);

float calcGilmore12_fixed(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file = 
      facilities::commonUtilities::joinPath(datadir, "opdep_fixed_Gilmore2012.dat");

   static AsciiTableModel gilmore_fixed(od_file);
   /// The Gilmore et al 2012 tables use MeV, so convert back from GeV.
   return gilmore_fixed.value(energy*1e3, redshift);
}

float calcGilmore12_fiducial(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file = 
      facilities::commonUtilities::joinPath(datadir, "opdep_fiducial_Gilmore2012.dat");
   static AsciiTableModel gilmore_fiducial(od_file);
   /// The Gilmore et al 2012 tables use MeV, so convert back from GeV.
   return gilmore_fiducial.value(energy*1e3, redshift);
}

float calcGeneric (float energy, float redshift){

//Parametric representation of tau(E,z) from Justin Finke


double E_0 = 80./(1.+redshift);    //energy in units of GeV
double E_1 = 450./(1.+redshift);
double E = (double)energy;

float tau;


if (E < 10. || redshift <= 0.) return 0.;

if (E < E_0)
   tau = (float)(2.6 * (1.+redshift) * pow(E_0/E_1, 2.) * pow(E/E_0, 4.));
   else if (energy < E_1)
        tau = (float)(2.6 * (1.+redshift) * pow(E/E_1, 2.));
		else
		   tau = (float)(2.6 * (1+redshift) * pow(E/E_1, 0.6));
		   
		 
		   
 return tau;
 
 }
		   

float calcSalamonStecker(float energy, float redshift){
// EBL model: Salamon & Stecker (ApJ 1998, 493:547-554)
//We are using here the model with metallicity correction (see paper)
//The paper has opacities up to z=3, for z>3 opacity remains constant according to Stecker



  int zindex=0,eindex=0;
  int i;
  float tau1,tau2;

  float zvalue[6] = {0., 0.1, 0.5, 1., 2., 3.};

  //evalues in units of GeV
  float evalue[14] = {10., 15., 20., 30., 40., 50., 60., 80., 100., 150., 200., 300., 400., 500.};

  /*From Fig.6 of the paper mentioned above */
  float tau[14][6] = {{0., 0., 0., 0., 0.023, 0.18},
		      {0., 0., 0., 0., 0.23, 0.8},
		      {0., 0., 0., 0.055, 0.8, 1.9},
		      {0., 0., 0.03, 0.31, 2.05, 3.8},
		      {0., 0., 0.11, 0.8, 3.4, 5.9},
		      {0., 0.01, 0.23, 1.23, 4.4, 7.},
		      {0., 0.02, 0.4, 1.7, 5.8, 8.},
		      {0., 0.05, 0.72, 2.8, 7.8, 10.0},
		      {0., 0.09, 1.1, 3.4, 9., 12.},
		      {0., 0.2, 2., 5.4, 10.7, 15.},
		      {0., 0.3, 2.6, 7., 13., 16.},
		      {0., 0.5, 3.9, 9., 15., 18.},
		      {0., 0.7, 4.9, 10., 17., 19.5},
		      {0., 0.8, 5.1, 10.1, 18.0, 20.}};


  if(redshift < 0.){
     std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
     redshift = 0.;
     }

//According to the model by the authors (Stecker and Salamon) opacities for z>3 are
//the same as for z=3. This is related to star formation among other things
//This is a really sensitive topic among theory people...  beware :)
  if (redshift > 3.) redshift = 3.;


  //Determine redshift index...
  if(redshift>=3.) zindex=5;
    else
      for(i=0;i<5;i++){
         if(redshift>=zvalue[i] && redshift<zvalue[i+1]){
	   zindex=i;
         }
       }

//Ebl attenuation is negligible for photons with E < 10 GeV (0.01 TeV)
//If energy > 500 GeV assume attenuation equivalent to energy = 500 GeV (outside Glast Energy Range anyway...)
  if(energy <= 10.) return 0;
    else if (energy >= 500.) {
      eindex = 13;
      energy = 500.;
      }
     else
       for(i=0;i<13;i++){
        if(energy>=evalue[i] && energy<evalue[i+1]){
	  eindex=i;
         }
        }


// Extrapolate in Energy for redshifts above and below the source
tau1 = tau[eindex][zindex]+(tau[eindex+1][zindex]-tau[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
tau2 = tau[eindex][zindex+1]+(tau[eindex+1][zindex+1]-tau[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);

//Extrapolate in redshift
return tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);

}



float calcGilmore(float energy, float redshift){


int zindex=0, eindex=-1;

float zvalue[16] = {0.0, 0.1, 0.5, 0.75, 1., 1.4, 1.8, 2., 2.4, 3., 4., 5., 6., 7., 8., 9.};

float evalue[61] = {1.000E+00, 1.122E+00, 1.259E+00, 1.413E+00, 1.585E+00, 1.778E+00, 1.995E+00, 2.239E+00, 2.512E+00, 2.818E+00, 3.162E+00, 
                    3.548E+00, 3.981E+00, 4.467E+00, 5.012E+00, 5.623E+00, 6.310E+00, 7.079E+00, 7.943E+00, 8.913E+00, 1.000E+01, 1.122E+01,
                    1.259E+01, 1.413E+01, 1.585E+01, 1.778E+01, 1.995E+01, 2.239E+01, 2.512E+01, 2.818E+01, 3.162E+01, 3.548E+01, 3.981E+01, 
					4.467E+01, 5.012E+01, 5.623E+01, 6.310E+01, 7.079E+01, 7.943E+01, 8.913E+01, 1.000E+02, 1.122E+02, 1.259E+02, 1.413E+02, 
					1.585E+02, 1.778E+02, 1.995E+02, 2.239E+02, 2.512E+02, 2.818E+02, 3.162E+02, 3.548E+02, 3.981E+02, 4.467E+02, 5.012E+02, 
					5.623E+02, 6.310E+02, 7.079E+02, 7.943E+02, 8.913E+02, 1.122E+03};

//Number of energy entries in the opacity table
  int MAXEINDEX = 61;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 16;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 9.){
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 9.0"<<std::endl;
#endif
       redshift=9.;
       }
if (energy >= 1.122e3) {
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for E < 1.122 TeV"<<std::endl;
#endif
       energy=1.122e3;
   } else if (energy < evalue[0]) return 0.;


float tautable[61][16]= 
{{0.0000E+00,1.9485E-06,1.4155E-05,2.5517E-05,3.9393E-05,6.4928E-05,9.0414E-05,1.0156E-04,1.1904E-04,1.3423E-04,1.4373E-04,1.4687E-04,1.4803E-04,1.4865E-04,1.4900E-04,1.4916E-04},
{0.0000E+00,2.0763E-06,1.5317E-05,2.7772E-05,4.3035E-05,7.1260E-05,9.9480E-05,1.1177E-04,1.3088E-04,1.4712E-04,1.5693E-04,1.6036E-04,1.6162E-04,1.6224E-04,1.6258E-04,1.6274E-04},
{0.0000E+00,2.2222E-06,1.6607E-05,3.0225E-05,4.6972E-05,7.8134E-05,1.0936E-04,1.2291E-04,1.4378E-04,1.6113E-04,1.7176E-04,1.7609E-04,1.7759E-04,1.7823E-04,1.7856E-04,1.7872E-04},
{0.0000E+00,2.3884E-06,1.8015E-05,3.2861E-05,5.1197E-05,8.5571E-05,1.2011E-04,1.3502E-04,1.5780E-04,1.7664E-04,1.8944E-04,1.9575E-04,1.9759E-04,1.9826E-04,1.9860E-04,1.9877E-04},
{0.0000E+00,2.5758E-06,1.9527E-05,3.5665E-05,5.5703E-05,9.3593E-05,1.3175E-04,1.4815E-04,1.7320E-04,1.9486E-04,2.1271E-04,2.2202E-04,2.2431E-04,2.2505E-04,2.2541E-04,2.2559E-04},
{0.0000E+00,2.7833E-06,2.1127E-05,3.8619E-05,6.0482E-05,1.0221E-04,1.4437E-04,1.6257E-04,1.9115E-04,2.1846E-04,2.4566E-04,2.5884E-04,2.6170E-04,2.6253E-04,2.6294E-04,2.6314E-04},
{0.0000E+00,3.0086E-06,2.2795E-05,4.1705E-05,6.5525E-05,1.1144E-04,1.5846E-04,1.7939E-04,2.1413E-04,2.5250E-04,2.9343E-04,3.1130E-04,3.1491E-04,3.1593E-04,3.1642E-04,3.1665E-04},
{0.0000E+00,3.2478E-06,2.4515E-05,4.4907E-05,7.0814E-05,1.2149E-04,1.7558E-04,2.0101E-04,2.4680E-04,3.0332E-04,3.6204E-04,3.8553E-04,3.9018E-04,3.9148E-04,3.9209E-04,3.9258E-04},
{0.0000E+00,3.4966E-06,2.6269E-05,4.8204E-05,7.6353E-05,1.3326E-04,1.9852E-04,2.3163E-04,2.9580E-04,3.7755E-04,4.5804E-04,4.8833E-04,4.9445E-04,4.9618E-04,4.9755E-04,5.0080E-04},
{0.0000E+00,3.7504E-06,2.8040E-05,5.1580E-05,8.2370E-05,1.4852E-04,2.3207E-04,2.7746E-04,3.6794E-04,4.8176E-04,5.8836E-04,6.2702E-04,6.3520E-04,6.3864E-04,6.4770E-04,6.6137E-04},
{0.0000E+00,4.0047E-06,2.9811E-05,5.5125E-05,8.9583E-05,1.7036E-04,2.8226E-04,3.4495E-04,4.6963E-04,6.2250E-04,7.6033E-04,8.0946E-04,8.2205E-04,8.4123E-04,8.8006E-04,9.1634E-04},
{0.0000E+00,4.2557E-06,3.1595E-05,5.9205E-05,9.9310E-05,2.0325E-04,3.5506E-04,4.4016E-04,6.0697E-04,8.0652E-04,9.8192E-04,1.0453E-03,1.0877E-03,1.1716E-03,1.2780E-03,1.3514E-03},
{0.0000E+00,4.4995E-06,3.3521E-05,6.4551E-05,1.1382E-04,2.5199E-04,4.5602E-04,5.6867E-04,7.8598E-04,1.0411E-03,1.2621E-03,1.3709E-03,1.5482E-03,1.7875E-03,2.0088E-03,2.1346E-03},
{0.0000E+00,4.7339E-06,3.5895E-05,7.2478E-05,1.3608E-04,3.2133E-04,5.9027E-04,7.3590E-04,1.0131E-03,1.3345E-03,1.6194E-03,1.9426E-03,2.4725E-03,2.9887E-03,3.3757E-03,3.5684E-03},
{0.0000E+00,4.9665E-06,3.9296E-05,8.4856E-05,1.6934E-04,4.1572E-04,7.6299E-04,9.4765E-04,1.2958E-03,1.6965E-03,2.1596E-03,3.1389E-03,4.3339E-03,5.2636E-03,5.8661E-03,6.1387E-03},
{0.0000E+00,5.2275E-06,4.4653E-05,1.0389E-04,2.1693E-04,5.3948E-04,9.7989E-04,1.2106E-03,1.6426E-03,2.1524E-03,3.2316E-03,5.6055E-03,7.8380E-03,9.3191E-03,1.0182E-02,1.0545E-02},
{0.0000E+00,5.5724E-06,5.3162E-05,1.3198E-04,2.8217E-04,6.9718E-04,1.2477E-03,1.5326E-03,2.0676E-03,2.8430E-03,5.5244E-03,1.0260E-02,1.3925E-02,1.6085E-02,1.7245E-02,1.7708E-02},
{0.0000E+00,6.1020E-06,6.6253E-05,1.7167E-04,3.6838E-04,8.9404E-04,1.5746E-03,1.9254E-03,2.6391E-03,4.1920E-03,1.0058E-02,1.8261E-02,2.3736E-02,2.6684E-02,2.8174E-02,2.8742E-02},
{0.0000E+00,6.9615E-06,8.5538E-05,2.2554E-04,4.7912E-04,1.1363E-03,1.9782E-03,2.4390E-03,3.6202E-03,6.9820E-03,1.8120E-02,3.0924E-02,3.8548E-02,4.2376E-02,4.4219E-02,4.4899E-02},
{0.0000E+00,8.3432E-06,1.1278E-04,2.9640E-04,6.1845E-04,1.4338E-03,2.5463E-03,3.2765E-03,5.5843E-03,1.2328E-02,3.1161E-02,4.9635E-02,5.9694E-02,6.4475E-02,6.6692E-02,6.7489E-02},
{0.0000E+00,1.0496E-05,1.4991E-04,3.8735E-04,7.9137E-04,1.8296E-03,3.5538E-03,4.9071E-03,9.4309E-03,2.1606E-02,5.0681E-02,7.5768E-02,8.8496E-02,9.4290E-02,9.6898E-02,9.7822E-02},
{0.0000E+00,1.3725E-05,1.9909E-04,5.0206E-04,1.0088E-03,2.4840E-03,5.5346E-03,8.0934E-03,1.6339E-02,3.6337E-02,7.8115E-02,1.1061E-01,1.2620E-01,1.3306E-01,1.3609E-01,1.3715E-01},
{0.0000E+00,1.8401E-05,2.6278E-04,6.4726E-04,1.3213E-03,3.7413E-03,9.2812E-03,1.3858E-02,2.7682E-02,5.8039E-02,1.1475E-01,1.5532E-01,1.7394E-01,1.8193E-01,1.8543E-01,1.8666E-01},
{0.0000E+00,2.4969E-05,3.4441E-04,8.5308E-04,1.8816E-03,6.1457E-03,1.5810E-02,2.3415E-02,4.4900E-02,8.8106E-02,1.6169E-01,2.1089E-01,2.3273E-01,2.4198E-01,2.4601E-01,2.4742E-01},
{0.0000E+00,3.3940E-05,4.5728E-04,1.2172E-03,2.9837E-03,1.0435E-02,2.6290E-02,3.8070E-02,6.9374E-02,1.2773E-01,2.1983E-01,2.7825E-01,3.0365E-01,3.1431E-01,3.1894E-01,3.2056E-01},
{0.0000E+00,4.5885E-05,6.5354E-04,1.9382E-03,5.0685E-03,1.7507E-02,4.1945E-02,5.9108E-02,1.0233E-01,1.7786E-01,2.8997E-01,3.5847E-01,3.8788E-01,4.0014E-01,4.0545E-01,4.0729E-01},
{0.0000E+00,6.2145E-05,1.0502E-03,3.3267E-03,8.7133E-03,2.8357E-02,6.3954E-02,8.7698E-02,1.4477E-01,2.3928E-01,3.7306E-01,4.5294E-01,4.8690E-01,5.0096E-01,5.0702E-01,5.0910E-01},
{0.0000E+00,9.0485E-05,1.8352E-03,5.7998E-03,1.4600E-02,4.3988E-02,9.3374E-02,1.2483E-01,1.9752E-01,3.1270E-01,4.7057E-01,5.6347E-01,6.0257E-01,6.1863E-01,6.2550E-01,6.2786E-01},
{0.0000E+00,1.5307E-04,3.2623E-03,9.8531E-03,2.3466E-02,6.5345E-02,1.3109E-01,1.7132E-01,2.6125E-01,3.9919E-01,5.8463E-01,6.9238E-01,7.3720E-01,7.5546E-01,7.6323E-01,7.6588E-01},
{0.0000E+00,2.8491E-04,5.6339E-03,1.6017E-02,3.6035E-02,9.3252E-02,1.7784E-01,2.2784E-01,3.3684E-01,5.0058E-01,7.1799E-01,8.4249E-01,8.9364E-01,9.1430E-01,9.2304E-01,9.2601E-01},
{0.0000E+00,5.2946E-04,9.2744E-03,2.4819E-02,5.2977E-02,1.2840E-01,2.3429E-01,2.9521E-01,4.2581E-01,6.1962E-01,8.7403E-01,1.0172E+00,1.0753E+00,1.0985E+00,1.1083E+00,1.1116E+00},
{0.0000E+00,9.3600E-04,1.4508E-02,3.6756E-02,7.4877E-02,1.7143E-01,3.0145E-01,3.7489E-01,5.3072E-01,7.5994E-01,1.0567E+00,1.2204E+00,1.2860E+00,1.3122E+00,1.3231E+00,1.3267E+00},
{0.0000E+00,1.5548E-03,2.1641E-02,5.2273E-02,1.0226E-01,2.2311E-01,3.8111E-01,4.6934E-01,6.5514E-01,9.2593E-01,1.2705E+00,1.4567E+00,1.5305E+00,1.5596E+00,1.5718E+00,1.5758E+00},
{0.0000E+00,2.4344E-03,3.0960E-02,7.1783E-02,1.3570E-01,2.8483E-01,4.7613E-01,5.8211E-01,8.0359E-01,1.1227E+00,1.5206E+00,1.7311E+00,1.8137E+00,1.8461E+00,1.8595E+00,1.8639E+00},
{0.0000E+00,3.6222E-03,4.2743E-02,9.5735E-02,1.7602E-01,3.5889E-01,5.9043E-01,7.1774E-01,9.8145E-01,1.3558E+00,1.8126E+00,2.0493E+00,2.1412E+00,2.1771E+00,2.1917E+00,2.1966E+00},
{0.0000E+00,5.1642E-03,5.7289E-02,1.2478E-01,2.2469E-01,4.4860E-01,7.2894E-01,8.8172E-01,1.1948E+00,1.6315E+00,2.1527E+00,2.4173E+00,2.5191E+00,2.5585E+00,2.5746E+00,2.5798E+00},
{0.0000E+00,7.1075E-03,7.5020E-02,1.6008E-01,2.8406E-01,5.5830E-01,8.9747E-01,1.0803E+00,1.4501E+00,1.9565E+00,2.5475E+00,2.8416E+00,2.9537E+00,2.9968E+00,3.0142E+00,3.0198E+00},
{0.0000E+00,9.5098E-03,9.6675E-02,2.0347E-01,3.5731E-01,6.9325E-01,1.1025E+00,1.3202E+00,1.7546E+00,2.3380E+00,3.0038E+00,3.3287E+00,3.4515E+00,3.4983E+00,3.5171E+00,3.5231E+00},
{0.0000E+00,1.2453E-02,1.2349E-01,2.5754E-01,4.4853E-01,8.5952E-01,1.3512E+00,1.6089E+00,2.1157E+00,2.7833E+00,3.5285E+00,3.8854E+00,4.0191E+00,4.0696E+00,4.0895E+00,4.0959E+00},
{0.0000E+00,1.6085E-02,1.5725E-01,3.2565E-01,5.6263E-01,1.0638E+00,1.6510E+00,1.9540E+00,2.5413E+00,3.2997E+00,4.1287E+00,4.5183E+00,4.6627E+00,4.7167E+00,4.7378E+00,4.7444E+00},
{0.0000E+00,2.0669E-02,2.0031E-01,4.1184E-01,7.0514E-01,1.3132E+00,2.0097E+00,2.3632E+00,3.0387E+00,3.8946E+00,4.8108E+00,5.2333E+00,5.3879E+00,5.4451E+00,5.4674E+00,5.4745E+00},
{0.0000E+00,2.6595E-02,2.5548E-01,5.2063E-01,8.8200E-01,1.6150E+00,2.4348E+00,2.8441E+00,3.6152E+00,4.5747E+00,5.5807E+00,6.0352E+00,6.1993E+00,6.2598E+00,6.2834E+00,6.2908E+00},
{0.0000E+00,3.4333E-02,3.2595E-01,6.5689E-01,1.0995E+00,1.9766E+00,2.9335E+00,3.4036E+00,4.2778E+00,5.3460E+00,6.4425E+00,6.9271E+00,7.1006E+00,7.1646E+00,7.1894E+00,7.1972E+00},
{0.0000E+00,4.4453E-02,4.1512E-01,8.2574E-01,1.3638E+00,2.4045E+00,3.5123E+00,4.0481E+00,5.0320E+00,6.2128E+00,7.3980E+00,7.9109E+00,8.0942E+00,8.1616E+00,8.1878E+00,8.1961E+00},
{0.0000E+00,5.7605E-02,5.2663E-01,1.0324E+00,1.6810E+00,2.9049E+00,4.1769E+00,4.7831E+00,5.8821E+00,7.1766E+00,8.4465E+00,8.9874E+00,9.1804E+00,9.2514E+00,9.2791E+00,9.2886E+00},
{0.0000E+00,7.4479E-02,6.6418E-01,1.2817E+00,2.0562E+00,3.4827E+00,4.9314E+00,5.6117E+00,6.8290E+00,8.2356E+00,9.5860E+00,1.0154E+01,1.0357E+01,1.0432E+01,1.0464E+01,1.0475E+01},
{0.0000E+00,9.5805E-02,8.3138E-01,1.5781E+00,2.4939E+00,4.1419E+00,5.7777E+00,6.5342E+00,7.8701E+00,9.3844E+00,1.0812E+01,1.1407E+01,1.1620E+01,1.1703E+01,1.1742E+01,1.1758E+01},
{0.0000E+00,1.2237E-01,1.0314E+00,1.9249E+00,2.9974E+00,4.8839E+00,6.7136E+00,7.5467E+00,8.9984E+00,1.0615E+01,1.2116E+01,1.2737E+01,1.2967E+01,1.3067E+01,1.3119E+01,1.3142E+01},
{0.0000E+00,1.5492E-01,1.2667E+00,2.3249E+00,3.5685E+00,5.7070E+00,7.7325E+00,8.6406E+00,1.0204E+01,1.1918E+01,1.3487E+01,1.4139E+01,1.4403E+01,1.4535E+01,1.4608E+01,1.4640E+01},
{0.0000E+00,1.9411E-01,1.5391E+00,2.7792E+00,4.2065E+00,6.6050E+00,8.8236E+00,9.8039E+00,1.1473E+01,1.3277E+01,1.4912E+01,1.5619E+01,1.5949E+01,1.6133E+01,1.6234E+01,1.6278E+01},
{0.0000E+00,2.4040E-01,1.8494E+00,3.2868E+00,4.9070E+00,7.5672E+00,9.9728E+00,1.1022E+01,1.2789E+01,1.4677E+01,1.6389E+01,1.7200E+01,1.7642E+01,1.7899E+01,1.8036E+01,1.8093E+01},
{0.0000E+00,2.9415E-01,2.1971E+00,3.8443E+00,5.6621E+00,8.5799E+00,1.1163E+01,1.2276E+01,1.4133E+01,1.6099E+01,1.7932E+01,1.8933E+01,1.9540E+01,1.9889E+01,2.0069E+01,2.0140E+01},
{0.0000E+00,3.5555E-01,2.5799E+00,4.4448E+00,6.4601E+00,9.6266E+00,1.2375E+01,1.3546E+01,1.5486E+01,1.7540E+01,1.9586E+01,2.0889E+01,2.1713E+01,2.2173E+01,2.2403E+01,2.2493E+01},
{0.0000E+00,4.2448E-01,2.9927E+00,5.0781E+00,7.2869E+00,1.0688E+01,1.3588E+01,1.4811E+01,1.6832E+01,1.9016E+01,2.1432E+01,2.3156E+01,2.4246E+01,2.4838E+01,2.5127E+01,2.5237E+01},
{0.0000E+00,5.0041E-01,3.4284E+00,5.7323E+00,8.1263E+00,1.1744E+01,1.4780E+01,1.6054E+01,1.8178E+01,2.0578E+01,2.3567E+01,2.5830E+01,2.7239E+01,2.7985E+01,2.8342E+01,2.8476E+01},
{0.0000E+00,5.8220E-01,3.8778E+00,6.3937E+00,8.9612E+00,1.2775E+01,1.5942E+01,1.7277E+01,1.9554E+01,2.2315E+01,2.6099E+01,2.9023E+01,3.0803E+01,3.1727E+01,3.2161E+01,3.2321E+01},
{0.0000E+00,6.6822E-01,4.3317E+00,7.0483E+00,9.7744E+00,1.3767E+01,1.7080E+01,1.8505E+01,2.1032E+01,2.4333E+01,2.9147E+01,3.2851E+01,3.5057E+01,3.6184E+01,3.6706E+01,3.6898E+01},
{0.0000E+00,7.5661E-01,4.7800E+00,7.6825E+00,1.0551E+01,1.4718E+01,1.8232E+01,1.9798E+01,2.2708E+01,2.6750E+01,3.2834E+01,3.7437E+01,4.0127E+01,4.1487E+01,4.2113E+01,4.2344E+01},
{0.0000E+00,8.4550E-01,5.2136E+00,8.2842E+00,1.1283E+01,1.5649E+01,1.9464E+01,2.1241E+01,2.4687E+01,2.9689E+01,3.7284E+01,4.2910E+01,4.6150E+01,4.7784E+01,4.8537E+01,4.8817E+01},
{0.0000E+00,9.3303E-01,5.6236E+00,8.8462E+00,1.1975E+01,1.6606E+01,2.0860E+01,2.2932E+01,2.7087E+01,3.3276E+01,4.2624E+01,4.9411E+01,5.3282E+01,5.5246E+01,5.6155E+01,5.6494E+01},
{0.0000E+00,9.8134E-01,6.0241E+00,9.6209E+00,1.3217E+01,1.8944E+01,2.5000E+01,2.8127E+01,3.5046E+01,4.6630E+01,6.4143E+01,7.5991E+01,8.2091E+01,1.0000E+10,1.0000E+10,1.0000E+10}};


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;


}




float calcFinke(float energy, float redshift){

//convert energy to TeV
energy /= 1.0e3;

int zindex=0, eindex=-1;

float zvalue[11]={0., 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 4.99};

float evalue[61]={1.000E-03, 1.122E-03, 1.259E-03, 1.413E-03, 1.585E-03, 1.778E-03, 1.995E-03, 2.239E-03, 2.512E-03, 2.818E-03, 3.162E-03,
                3.548E-03, 3.981E-03, 4.467E-03, 5.012E-03, 5.623E-03, 6.310E-03, 7.079E-03, 7.943E-03, 8.913E-03, 1.000E-02, 1.122E-02, 
                1.259E-02, 1.413E-02, 1.585E-02, 1.778E-02, 1.995E-02, 2.239E-02, 2.512E-02, 2.818E-02, 3.162E-02, 3.548E-02, 3.981E-02,
				4.467E-02, 5.012E-02, 5.623E-02, 6.310E-02, 7.079E-02, 7.943E-02, 8.913E-02, 1.000E-01, 1.122E-01, 1.259E-01, 1.413E-01, 
                1.585E-01, 1.778E-01, 1.995E-01, 2.239E-01, 2.512E-01, 2.818E-01, 3.162E-01, 3.548E-01, 3.981E-01, 4.467E-01, 5.012E-01, 
				5.623E-01, 6.310E-01, 7.079E-01, 7.943E-01, 8.913E-01, 1.000E+00};

//Number of energy entries in the opacity table
  int MAXEINDEX = 61;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 11;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 4.99){
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 4.99"<<std::endl;
#endif
       redshift=4.99;
       }
if (energy >= 1.) {
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for E < 1. TeV"<<std::endl;
#endif
       energy=1.;
   } else if (energy < evalue[0]) return 0.;


float tautable [61][11] = 
{{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.08E-59},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,2.29E-57},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,6.74E-58,2.19E-56},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.68E-56,3.36E-07},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,9.78E-58,1.84E-09,1.29E-05},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.71E-56,4.79E-06,1.15E-04},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,3.66E-59,1.01E-08,7.71E-05,5.06E-04},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,6.32E-57,1.03E-05,4.77E-04,1.51E-03},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,8.72E-58,5.68E-56,1.47E-04,1.69E-03,3.65E-03},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.78E-56,1.18E-06,8.49E-04,4.42E-03,7.67E-03},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,1.00E-100,1.16E-57,1.86E-09,4.42E-05,2.82E-03,9.62E-03,1.45E-02},
{0.000E+00,1.00E-100,1.00E-100,1.00E-100,5.11E-58,2.01E-56,6.06E-06,3.84E-04,7.01E-03,1.85E-02,2.55E-02},
{0.000E+00,1.00E-100,1.00E-100,3.03E-59,1.36E-56,4.69E-09,9.54E-05,1.61E-03,1.46E-02,3.26E-02,4.20E-02},
{0.000E+00,1.00E-100,1.00E-100,5.52E-57,2.19E-12,4.96E-06,5.85E-04,4.52E-03,2.72E-02,5.34E-02,6.53E-02},
{0.000E+00,1.00E-100,2.44E-57,4.86E-56,1.92E-06,7.14E-05,2.01E-03,1.02E-02,4.62E-02,8.24E-02,9.73E-02},
{0.000E+00,3.10E-58,2.76E-56,2.29E-07,3.49E-05,4.10E-04,5.11E-03,1.99E-02,7.35E-02,1.22E-01,1.40E-01},
{0.000E+00,7.80E-57,2.89E-08,9.05E-06,2.32E-04,1.38E-03,1.09E-02,3.52E-02,1.11E-01,1.73E-01,1.95E-01},
{0.000E+00,4.52E-56,3.00E-06,8.05E-05,8.27E-04,3.46E-03,2.04E-02,5.76E-02,1.59E-01,2.40E-01,2.65E-01},
{0.000E+00,3.84E-07,3.24E-05,3.46E-04,2.16E-03,7.29E-03,3.50E-02,8.85E-02,2.22E-01,3.23E-01,3.52E-01},
{0.000E+00,6.47E-06,1.57E-04,9.88E-04,4.66E-03,1.37E-02,5.59E-02,1.30E-01,3.02E-01,4.27E-01,4.60E-01},
{0.000E+00,3.92E-05,4.74E-04,2.27E-03,8.89E-03,2.35E-02,8.45E-02,1.83E-01,4.00E-01,5.55E-01,5.92E-01},
{0.000E+00,1.25E-04,1.12E-03,4.50E-03,1.55E-02,3.77E-02,1.22E-01,2.51E-01,5.22E-01,7.12E-01,7.51E-01},
{0.000E+00,2.96E-04,2.25E-03,8.06E-03,2.52E-02,5.71E-02,1.71E-01,3.36E-01,6.72E-01,9.01E-01,9.44E-01},
{0.000E+00,5.97E-04,4.07E-03,1.34E-02,3.85E-02,8.29E-02,2.33E-01,4.41E-01,8.53E-01,1.13E+00,1.17E+00},
{0.000E+00,1.09E-03,6.82E-03,2.09E-02,5.64E-02,1.16E-01,3.10E-01,5.72E-01,1.07E+00,1.41E+00,1.45E+00},
{0.000E+00,1.83E-03,1.07E-02,3.11E-02,7.99E-02,1.59E-01,4.07E-01,7.32E-01,1.34E+00,1.74E+00,1.78E+00},
{0.000E+00,2.90E-03,1.61E-02,4.46E-02,1.10E-01,2.13E-01,5.27E-01,9.28E-01,1.65E+00,2.13E+00,2.17E+00},
{0.000E+00,4.36E-03,2.32E-02,6.21E-02,1.49E-01,2.82E-01,6.76E-01,1.17E+00,2.03E+00,2.60E+00,2.64E+00},
{0.000E+00,6.32E-03,3.25E-02,8.47E-02,1.98E-01,3.68E-01,8.58E-01,1.45E+00,2.48E+00,3.16E+00,3.18E+00},
{0.000E+00,8.91E-03,4.46E-02,1.14E-01,2.60E-01,4.75E-01,1.08E+00,1.80E+00,3.02E+00,3.81E+00,3.82E+00},
{0.000E+00,1.23E-02,6.02E-02,1.51E-01,3.38E-01,6.09E-01,1.35E+00,2.22E+00,3.64E+00,4.58E+00,4.56E+00},
{0.000E+00,1.67E-02,8.04E-02,1.98E-01,4.36E-01,7.75E-01,1.69E+00,2.72E+00,4.38E+00,5.45E+00,5.42E+00},
{0.000E+00,2.25E-02,1.06E-01,2.59E-01,5.60E-01,9.80E-01,2.08E+00,3.31E+00,5.23E+00,6.46E+00,6.39E+00},
{0.000E+00,3.00E-02,1.40E-01,3.35E-01,7.13E-01,1.23E+00,2.56E+00,4.00E+00,6.20E+00,7.59E+00,7.49E+00},
{0.000E+00,3.97E-02,1.83E-01,4.31E-01,9.03E-01,1.54E+00,3.12E+00,4.80E+00,7.29E+00,8.84E+00,8.71E+00},
{0.000E+00,5.23E-02,2.37E-01,5.51E-01,1.13E+00,1.90E+00,3.78E+00,5.71E+00,8.51E+00,1.02E+01,1.00E+01},
{0.000E+00,6.84E-02,3.05E-01,6.99E-01,1.41E+00,2.34E+00,4.54E+00,6.74E+00,9.85E+00,1.17E+01,1.15E+01},
{0.000E+00,8.90E-02,3.90E-01,8.79E-01,1.75E+00,2.84E+00,5.39E+00,7.89E+00,1.13E+01,1.33E+01,1.30E+01},
{0.000E+00,1.15E-01,4.94E-01,1.10E+00,2.14E+00,3.42E+00,6.35E+00,9.13E+00,1.28E+01,1.50E+01,1.47E+01},
{0.000E+00,1.46E-01,6.19E-01,1.35E+00,2.59E+00,4.08E+00,7.39E+00,1.05E+01,1.44E+01,1.67E+01,1.64E+01},
{0.000E+00,1.85E-01,7.67E-01,1.65E+00,3.09E+00,4.81E+00,8.52E+00,1.19E+01,1.61E+01,1.85E+01,1.81E+01},
{0.000E+00,2.31E-01,9.40E-01,1.98E+00,3.66E+00,5.61E+00,9.72E+00,1.33E+01,1.78E+01,2.04E+01,1.99E+01},
{0.000E+00,2.84E-01,1.14E+00,2.36E+00,4.28E+00,6.47E+00,1.10E+01,1.48E+01,1.95E+01,2.22E+01,2.18E+01},
{0.000E+00,3.45E-01,1.36E+00,2.76E+00,4.94E+00,7.37E+00,1.22E+01,1.63E+01,2.13E+01,2.41E+01,2.37E+01},
{0.000E+00,4.14E-01,1.60E+00,3.21E+00,5.63E+00,8.29E+00,1.35E+01,1.78E+01,2.30E+01,2.61E+01,2.57E+01},
{0.000E+00,4.89E-01,1.85E+00,3.67E+00,6.35E+00,9.23E+00,1.48E+01,1.93E+01,2.48E+01,2.81E+01,2.78E+01},
{0.000E+00,5.70E-01,2.13E+00,4.15E+00,7.08E+00,1.02E+01,1.60E+01,2.08E+01,2.66E+01,3.02E+01,2.99E+01},
{0.000E+00,6.56E-01,2.41E+00,4.64E+00,7.80E+00,1.11E+01,1.72E+01,2.23E+01,2.86E+01,3.24E+01,3.23E+01},
{0.000E+00,7.45E-01,2.70E+00,5.12E+00,8.51E+00,1.20E+01,1.84E+01,2.38E+01,3.07E+01,3.48E+01,3.48E+01},
{0.000E+00,8.36E-01,2.98E+00,5.60E+00,9.19E+00,1.29E+01,1.96E+01,2.53E+01,3.30E+01,3.74E+01,3.75E+01},
{0.000E+00,9.26E-01,3.26E+00,6.06E+00,9.85E+00,1.37E+01,2.09E+01,2.70E+01,3.55E+01,4.03E+01,4.06E+01},
{0.000E+00,1.02E+00,3.53E+00,6.50E+00,1.05E+01,1.45E+01,2.22E+01,2.89E+01,3.83E+01,4.35E+01,4.39E+01},
{0.000E+00,1.10E+00,3.79E+00,6.92E+00,1.11E+01,1.54E+01,2.36E+01,3.10E+01,4.15E+01,4.70E+01,4.76E+01}}; 
						   


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;


}

float calcFranceschini(float energy, float redshift){

//convert energy to TeV
energy /= 1.0e3;


int zindex=0, eindex=-1;

float zvalue[10] = {0., 0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0};

float evalue[26] = {0.01, 0.0200, 0.0240, 0.0289, 0.0347, 0.0417, 0.0502, 0.0603, 0.0726, 0.0873, 0.1040, 0.1260,
                    0.1510, 0.1820, 0.2190, 0.2630, 0.3160, 0.3810, 0.4580, 0.5500, 0.6620, 0.7960, 0.9570, 
					1.1500, 1.3800, 1.6600};

//Number of energy entries in the opacity table
  int MAXEINDEX = 26;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 10;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 3.){
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
#endif
       redshift=3.;
       }
if (energy >= 1.66) {
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
#endif
       energy=1.66;
   } else if (energy < evalue[0]) return 0.;

float tautable[26][10]={ {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                         {0., 0,0,0,0,0.0000021,0.004933,0.0399,0.1157,0.2596},
						 {0., 0,0,0,0,0.000188,0.01284,0.0718,0.1783,0.3635},
						 {0., 0,0,0,0,0.001304,0.0279,0.1188,0.2598,0.4919}, 
                         {0., 0,0,0,0.000488,0.004558,0.0533,0.1833,0.3635,0.6517}, 
                         {0., 0,0,5.25E-05,0.002276,0.01157,0.0921,0.2689,0.4967,0.8548}, 
						 {0., 0,9.45E-05,5.41E-04,0.006575,0.02436,0.148,0.3836,0.6745,1.118}, 
                         {0., 1.10E-04,4.24E-04,1.92E-03,0.014592,0.04512,0.2275,0.5434,0.9179,1.465}, 
                         {0., 3.09E-04,1.10E-03,4.55E-03,0.02771,0.07684,0.343,0.7707,1.251,1.917}, 
						 {0., 6.56E-04,2.26E-03,8.90E-03,0.04808,0.1248,0.5137,1.092,1.703,2.503}, 
                         {0., 1.21E-03,4.10E-03,1.58E-02,0.07958,0.1984,0.764,1.537,2.302,3.249}, 
						 {0., 2.11E-03,7.04E-03,2.69E-02,0.1284,0.3109,1.12,2.133,3.073,4.181}, 
                         {0., 3.53E-03,1.17E-02,4.41E-02,0.2031,0.478,1.607,2.905,4.042,5.318}, 
                         {0., 5.71E-03,1.87E-02,7.01E-02,0.3134,0.7163,2.247,3.875,5.225,6.673}, 
						 {0., 8.92E-03,2.91E-02,0.1082,0.4696,1.04,3.056,5.055,6.627,8.241}, 
						 {0., 1.35E-02,4.38E-02,0.1618,0.6809,1.461,4.042,6.438,8.226,9.997}, 
						 {0., 1.98E-02,6.37E-02,0.2338,0.9517,1.981,5.192,7.989,9.977,11.89}, 
                         {0., 2.79E-02,8.94E-02,0.3256,1.281,2.594,6.474,9.65,11.81,13.89}, 
                         {0., 3.80E-02,0.1205,0.4356,1.661,3.284,7.836,11.34,13.67,15.93}, 
                         {0., 4.96E-02,0.1563,0.5607,2.082,4.023,9.214,13.01,15.51,18.08}, 
                         {0., 6.23E-02,0.1953,0.6961,2.524,4.779,10.55,14.63,17.39,20.45}, 
                         {0., 7.58E-02,0.2364,0.8373,2.967,5.517,11.82,16.25,19.49,23.27}, 
						 {0., 8.92E-02,0.2768,0.975,3.389,6.21,13.03,18.04,22.02,26.81}, 
                         {0., 0.1019,0.3152,1.105,3.779,6.846,14.29,20.21,25.22,31.33}, 
                         {0., 0.1136,0.3501,1.223,4.129,7.432,15.73,22.98,29.37,37.23}, 
                         {0., 0.124,0.381,1.327,4.444,8.01,17.54,26.58,34.78,45.09}};

  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;

  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;


}

float calcKneiske(float energy, float redshift){
/************************************************************************
EBL model 0: Kneiske, Bretz, Mannheim, Hartmann (A&A 413, 807-815, 2004)
  Valid for redshift <= 5.0
  Here we have implemented the "best fit" model from their paper
************************************************************************/



int zindex=0, eindex=-1;


float zvalue[51];

for(int i=0; i<=50; i++) zvalue[i]= 0.1*i;

float evalue[11];

for(int i=0; i<=10; i++) evalue[i]= pow(10., 0.8+0.18*i);

//Number of energy entries in the opacity table
  int MAXEINDEX = 11;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 51;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 5.){
#ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
#endif
       redshift=5.;
       }
if (energy >= 350.) {
   return calcKneiske_extendedEnergyRange(energy, redshift);
   } else if (energy < evalue[0]) return 0.;


float tautable[11][51] = { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.002,0.003,0.005,0.007,0.009,0.011,0.014,0.018,0.02,0.021,0.026,0.028,0.03,0.033,0.036,0.038,0.041,0.044,0.044},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.003,0.007,0.012,0.018,0.026,0.031,0.041,0.047,0.056,0.064,0.075,0.082,0.093,0.101,0.11,0.122,0.131,0.137,0.142,0.156,0.161,0.17,0.176,0.182,0.188,0.197,0.198,0.2,0.209,0.217,0.221,0.219,0.224},
{0.,0.,0.,0.,0.,0.,0.,0.,0.002,0.007,0.016,0.031,0.047,0.063,0.082,0.105,0.122,0.148,0.165,0.195,0.215,0.234,0.263,0.28,0.301,0.326,0.346,0.362,0.381,0.401,0.423,0.437,0.446,0.468,0.489,0.508,0.515,0.53,0.543,0.557,0.573,0.58,0.599,0.609,0.613,0.616,0.619,0.64,0.653,0.66,0.662},
{0.,0.,0.,0.002,0.006,0.016,0.032,0.057,0.091,0.131,0.178,0.229,0.278,0.325,0.375,0.416,0.463,0.509,0.558,0.598,0.646,0.687,0.732,0.778,0.814,0.854,0.89,0.93,0.97,0.999,1.036,1.067,1.089,1.119,1.16,1.18,1.213,1.234,1.257,1.278,1.294,1.311,1.333,1.34,1.36,1.377,1.398,1.415,1.418,1.426,1.432},
{0.,0.004,0.014,0.032,0.063,0.104,0.158,0.222,0.297,0.382,0.475,0.572,0.67,0.763,0.854,0.946,1.036,1.125,1.212,1.296,1.377,1.46,1.54,1.611,1.687,1.755,1.828,1.89,1.958,2.021,2.069,2.126,2.178,2.231,2.286,2.325,2.363,2.4,2.436,2.475,2.503,2.539,2.569,2.587,2.602,2.618,2.658,2.676,2.698,2.71,2.725},
{0.,0.021,0.055,0.107,0.178,0.267,0.377,0.509,0.661,0.826,1.002,1.184,1.365,1.545,1.721,1.894,2.063,2.229,2.388,2.542,2.693,2.838,2.975,3.107,3.229,3.348,3.466,3.576,3.676,3.772,3.863,3.948,4.019,4.107,4.19,4.266,4.325,4.376,4.439,4.491,4.542,4.59,4.638,4.688,4.721,4.751,4.786,4.817,4.841,4.865,4.889},
{0.,0.053,0.132,0.244,0.394,0.584,0.815,1.086,1.388,1.711,2.048,2.39,2.729,3.053,3.371,3.666,3.958,4.242,4.513,4.775,5.02,5.26,5.492,5.709,5.916,6.113,6.308,6.483,6.652,6.811,6.969,7.119,7.25,7.379,7.505,7.618,7.72,7.822,7.921,8.01,8.098,8.179,8.251,8.307,8.363,8.414,8.468,8.517,8.564,8.599,8.624},
{0.,0.124,0.304,0.554,0.878,1.273,1.731,2.242,2.797,3.379,3.973,4.569,5.153,5.716,6.263,6.776,7.275,7.753,8.208,8.641,9.054,9.448,9.819,10.173,10.507,10.822,11.116,11.394,11.655,11.897,12.128,12.349,12.541,12.73,12.908,13.072,13.218,13.351,13.48,13.601,13.716,13.829,13.932,14.005,14.073,14.136,14.203,14.27,14.328,14.373,14.41},
{0.,0.287,0.686,1.207,1.855,2.616,3.469,4.399,5.38,6.387,7.397,8.391,9.346,10.253,11.114,11.91,12.671,13.384,14.054,14.686,15.273,15.828,16.346,16.829,17.289,17.71,18.106,18.472,18.817,19.142,19.445,19.732,19.979,20.22,20.459,20.676,20.873,21.055,21.238,21.407,21.565,21.726,21.876,21.993,22.095,22.191,22.298,22.41,22.508,22.59,22.657},
{0.,0.61,1.407,2.396,3.568,4.886,6.311,7.798,9.313,10.818,12.284,13.689,15.012,16.247,17.404,18.463,19.467,20.394,21.26,22.069,22.824,23.533,24.203,24.832,25.42,25.968,26.491,26.988,27.46,27.9,28.315,28.705,29.075,29.446,29.812,30.152,30.463,30.746,31.031,31.311,31.582,31.852,32.103,32.325,32.528,32.719,32.912,33.1,33.271,33.411,33.539}};


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;

}

float calcPrimack05(float energy, float redshift) {
   return Primack05::instance().value(energy, redshift);
}

// float calcPrimack05(float energy, float redshift){
// // EBL model 1: Primack & Bullock (2005)
   
//    int zindex=0, eindex=-1;
//    float tau1,tau2, **tautables, tauvalue;
   
//    float zvalue[17] = {0., 0.1, 0.25, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 
//                        4., 4.5, 5., 5.5, 6., 6.5, 7.};
   
//    float evalue[40]={2., 2.488, 3.0950, 3.8510, 4.7910, 5.9600, 7.4150, 
//                      9.2250, 1.1480e+1, 1.4280e+1, 1.7760e+1, 2.2100e+1, 
//                      2.7490e+1, 3.4200e+1, 4.2550e+1, 5.2930e+1, 6.5850e+1, 
//                      8.1920e+1, 1.0190e+2, 1.2680e+2, 1.5770e+2, 1.9620e+2,
//                      2.4410e+2, 3.0370e+2, 3.78e+2, 4.70e+2, 5.85e+2, 
//                      7.28e+2, 9.05e+2, 1.13e+3, 1.40e+3, 1.74e+3, 2.17e+3,
//                      2.70e+3, 3.36e+3, 4.18e+3, 5.19e+3, 6.46e+3, 8.04e+3,
//                      1.00e+4};

//   float tau0[1] ={0.};

//   float tau01[40] = {0., 0.0000001, 0.0000003, 0.0000008, 0.0000019, 0.0000042, 0.0000099, 0.0000238, 0.0000552, 0.0001199, 0.0002443, 0.0004892, 0.0010674, 0.0024132, 0.0050824, 0.0096015, 0.0163697, 0.0256939, 0.0380436, 0.0546355, 0.0777303, 0.1105819, 0.1570698, 0.2208418, 0.3037539, 0.4054053, 0.5221997, 0.6465960, 0.7691562, 0.8799724, 0.9709822, 1.0396483, 1.0923208, 1.1416358, 1.2064721, 1.3020910, 1.4366617, 1.6201703, 1.8610741, 2.1794928};

//   float tau025[40]={0.0000001, 0.0000005, 0.0000014, 0.0000035, 0.0000076, 0.0000171, 0.0000407, 0.0000968, 0.0002186, 0.0004617, 0.0009238, 0.0019085, 0.0042851, 0.0094652, 0.0190260, 0.0342809, 0.0560866, 0.0850821, 0.1230088, 0.1742253, 0.2459901, 0.3477387, 0.4895209, 0.6785036, 0.9165204, 1.1991670, 1.5127269, 1.8345708, 2.1389193, 2.4018006, 2.6075906, 2.7604884, 2.8892735, 3.0378014, 3.2507707, 3.5584134, 3.9840049, 4.5473202, 5.2793418, 6.2633874};

//   float tau05[40] = {0.0000010, 0.0000029, 0.0000067, 0.0000144, 0.0000312, 0.0000731, 0.0001749, 0.0004026, 0.0008668, 0.0017518, 0.0035597, 0.0078349, 0.0174737, 0.0360969, 0.0669131, 0.1117945, 0.1714340, 0.2479820, 0.3487295, 0.4873202, 0.6822443, 0.9530774, 1.3148751, 1.7737483, 2.3246879, 2.9468168, 3.6000372, 4.2326759, 4.7939800, 5.2462781, 5.5866135, 5.8646962, 6.1662335, 6.5838559, 7.1977531, 8.0583552, 9.1996407, 10.6771777, 12.6293777, 15.4169639};

//   float tau10[35] = {0.0000092, 0.0000193, 0.0000395, 0.0000864, 0.0002063, 0.0004962, 0.0011276, 0.0023725, 0.0047253, 0.0097697, 0.0218142, 0.0480538, 0.0964923, 0.1726325, 0.2780015, 0.4126927, 0.5828524, 0.8065691, 1.1128001, 1.5344083, 2.1009211, 2.8323235, 3.7352715, 4.7970329, 5.9754211, 7.1966803, 8.3668533, 9.3951736, 10.2286613, 10.8944157, 11.5071431, 12.2380499, 13.2666033, 14.7188514, 16.6776653};

//   float tau15[30] = {0.0000289, 0.0000579, 0.0001275, 0.0003103, 0.0007540, 0.0017071, 0.0035602, 0.0070077, 0.0144907, 0.0327386, 0.0725108, 0.1446813, 0.2558033, 0.4069964, 0.5984519, 0.8400723, 1.1570757, 1.5847601, 2.1601107, 2.9160698, 3.8756493, 5.0484552, 6.4230870, 7.9519874, 9.5477112, 11.0982863, 12.4918875, 13.6635642, 14.6535502, 15.6091534};

//   float tau20[27] = {0.0000625, 0.0001368, 0.0003343, 0.0008154, 0.0018518, 0.0038687, 0.0076352, 0.0157515, 0.0356483, 0.0792690, 0.1589536, 0.2827625, 0.4531616, 0.6713563, 0.9477445, 1.3074161, 1.7843432, 2.4149232, 3.2333568, 4.2665894,  5.5296426, 7.0184022, 8.6950112, 10.4807643, 12.2666076, 13.9356198, 15.4071088};

//   float tau25[26] = {0.0001093, 0.0002630, 0.0006461, 0.0014943, 0.0031830, 0.0063392, 0.0128773, 0.0287818, 0.0647205,
// 0.1332950, 0.2447116, 0.4043603, 0.6147429, 0.8837033, 1.2308171, 1.6848746, 2.2793325, 3.0487025, 4.0236373, 5.2262532, 6.6640363, 8.3175392, 10.1311691, 12.0153741, 13.8593746, 15.5607761};

//   float tau30[25] = {0.0002166, 0.0005354, 0.0012768, 0.0028085, 0.0057098, 0.0112874, 0.0242226, 0.0548200, 0.1172134, 0.2248061, 0.3857926, 0.6030031, 0.8803797, 1.2328916, 1.6865156, 2.2734022, 3.0271987, 3.9793934, 5.1548374, 6.5690724, 8.2179971, 10.0657856, 12.0419091, 14.0464631, 15.9680728};

//   float tau35[23] = {0.0005791, 0.0013993, 0.0031327, 0.0064560, 0.0127024, 0.0264943, 0.0592172, 0.1280992, 0.2507395, 0.4384755, 0.6944182, 1.0198036, 1.4263889, 1.9406097, 2.5965219, 3.4284109, 4.4663397, 5.7331447, 7.2439184, 8.9999333, 10.9746883, 13.1074494, 15.3066731};

//   float tau40[22] = {0.0015020, 0.0033793, 0.0070182, 0.0138902, 0.0289108, 0.0638794, 0.1374250, 0.2702333, 0.4775225, 0.7651981, 1.1348334, 1.5962845, 2.1734261, 2.9008888, 3.8146649, 4.9450719, 6.3128560, 7.9304825, 9.7980313, 11.8928097, 14.1623445, 16.5244052};

//   float tau45[20] = {0.0034667, 0.0072649, 0.0144116, 0.0297863, 0.0657044, 0.1418342, 0.2810762, 0.5021846, 0.8146987, 1.2220693, 1.7341971, 2.3738843, 3.1743106, 4.1729259, 5.4022331, 6.8829534, 8.6250875, 10.6260443, 12.8610481, 15.2794450};

//   float tau50[19] = {0.0071092, 0.0141624, 0.0288313, 0.0632518, 0.1382215, 0.2787677, 0.5066919, 0.8348414, 1.2687125, 1.8185965, 2.5079854, 3.3704304, 4.4431075, 5.7607600, 7.3463606, 9.2094083, 11.3455301, 13.7276341, 16.3016469};

//   float tau55[18] = {0.0132880, 0.0265899, 0.0574533, 0.1270139, 0.2628803, 0.4906496, 0.8268285, 1.2789404, 1.8567263, 2.5842328, 3.4973174, 4.6354043, 6.0346728, 7.7216602, 9.7080199, 11.9896983, 14.5392316, 17.2989585};

//   float tau60[16] = {0.0235960, 0.0496899, 0.1103218, 0.2348820, 0.4537640, 0.7883757, 1.2486188, 1.8424744, 2.5923866, 3.5366242, 4.7194613, 6.1811159, 7.9511663, 10.0444068, 12.4598110, 15.1731650};

//   float tau65[15] = {0.0414492, 0.0911960, 0.1991044, 0.3999018, 0.7218841, 1.1796636, 1.7796795, 2.5393540, 3.4965811, 4.7009047, 6.2001010, 8.0305576, 10.2107500, 12.7428649, 15.6087370};

//   float tau70[14] = {0.0724103, 0.1605235, 0.3354485, 0.6327967, 1.0749973, 1.6703149, 2.4294342, 3.3849692, 4.5888035,  6.0974870, 7.9579857, 10.1977161, 12.8222616, 15.8189200};

//   tautables = (float **)malloc(17*sizeof(float*));
//   tautables[0] = tau0;
//   tautables[1] = tau01;
//   tautables[2] = tau025;
//   tautables[3] = tau05;
//   tautables[4] = tau10;
//   tautables[5] = tau15;
//   tautables[6] = tau20;
//   tautables[7] = tau25;
//   tautables[8] = tau30;
//   tautables[9] = tau35;
//   tautables[10] = tau40;
//   tautables[11] = tau45;
//   tautables[12] = tau50;
//   tautables[13] = tau55;
//   tautables[14] = tau60;
//   tautables[15] = tau65;
//   tautables[16] = tau70;


//   int MaxEindex[17] = {0, 39, 39, 39, 34, 29, 26, 25, 24, 22, 21, 19, 18, 
//                        17, 15, 14, 13};

//   if (redshift < 0.) {
//      std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
//      redshift = 0.;
//   } else if (redshift > 7.) {
// #ifdef DEBUG
//      std::cerr<<"Maximum redshift for this model is z = 7.0... Calculating opacity for z = 7.0"<<std::endl;
// #endif
//      redshift=7.;
//   }

//   if (energy >= 10000.) {
//      std::cerr<<"This EBL model is only valid for E <= 10 TeV..."<<std::endl;
//      energy = 1e4;
//   } else if (energy < evalue[0]) {
//      return 0.;
//   }

//   //Determine redshift index...
//   for (int i=0; i<16; i++) {
//      if (redshift >= zvalue[i] && redshift < zvalue[i+1]) {
//         zindex = i;
//      }
//   }
//   if (redshift >= zvalue[16]) {
//      zindex = 16;
//   }

//   // Determine energy index
//   for (int i=0; i<39; i++) {
//      if (energy >= evalue[i] && energy < evalue[i+1]) {
//         eindex = i;
//      }
//   }
//   if (eindex < 0) {
//      return 0;
//   }

//   if (zindex < 16) { //Find tau for redshifts above and below source
//      if (eindex < MaxEindex[zindex]) {
//         tau1 = tautables[zindex][eindex]+(tautables[zindex][eindex+1]-tautables[zindex][eindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
//      } else {
//         tau1 = tautables[zindex][MaxEindex[zindex]];
//      }

//      if (eindex < MaxEindex[zindex+1]) {
//         tau2 = tautables[zindex+1][eindex]+(tautables[zindex+1][eindex+1]-tautables[zindex+1][eindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
//      } else {
//         tau2 = tautables[zindex+1][MaxEindex[zindex+1]];
//      }

//   //interpolate in redshift
//      tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
//   } else {
//     //Use tau for source at z = 7
//      if (eindex < MaxEindex[zindex]) {
//         tau1 = tautables[zindex][eindex]+(tautables[zindex][eindex+1]-tautables[zindex][eindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
//      } else {
//         tau1 = tautables[zindex][MaxEindex[zindex]];
//      }
//      tauvalue = tau1;
//   }

//   free(tautables);
//   return tauvalue;
// }

float calcKneiske_HighUV(float energy, float redshift){
/************************************************************************
EBL model 2: Kneiske, Bretz, Mannheim, Hartmann (A&A 413, 807-815, 2004)
  Valid for redshift <= 5.0
  Here we have implemented the "High UV" model from their paper
************************************************************************/



int zindex=0, eindex=-1;

float zvalue[51];

for(int i=0; i<=50; i++) zvalue[i]= 0.1*i;

float evalue[14];

for(int i=0; i<=13; i++) evalue[i]= pow(10., 0.26+0.18*i);

//Number of energy entries in the opacity table
  int MAXEINDEX = 14;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 51;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 5.){
       #ifdef DEBUG
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
       #endif
       redshift=5.;
       }
if (energy >= 350.) {
   return calcKneiske_HighUV_extendedEnergyRange(energy, redshift);
} else if (energy < evalue[0]) return 0.;

float tautable[14][51] = {
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.002,0.002},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.003,0.003,0.003,0.004,0.004,0.004,0.005,0.005,0.006,0.006,0.007,0.007,0.008,0.008,0.009,0.009,0.01,0.01,0.011,0.011,0.012,0.013,0.013,0.014},
{0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.003,0.004,0.005,0.006,0.007,0.009,0.01,0.011,0.013,0.014,0.016,0.017,0.019,0.021,0.023,0.025,0.027,0.029,0.031,0.033,0.036,0.038,0.04,0.043,0.045,0.047,0.05,0.052,0.054,0.056,0.058,0.061,0.063,0.064,0.066},
{0.,0.,0.,0.,0.,0.001,0.002,0.003,0.004,0.006,0.008,0.009,0.011,0.014,0.017,0.021,0.025,0.03,0.034,0.039,0.045,0.051,0.056,0.063,0.069,0.076,0.083,0.089,0.097,0.103,0.11,0.117,0.125,0.132,0.138,0.145,0.152,0.158,0.165,0.171,0.177,0.182,0.187,0.192,0.197,0.201,0.206,0.21,0.214,0.216,0.219},
{0.,0.,0.,0.001,0.002,0.004,0.006,0.01,0.014,0.021,0.031,0.043,0.056,0.071,0.086,0.102,0.117,0.133,0.15,0.167,0.184,0.201,0.218,0.234,0.251,0.267,0.283,0.298,0.314,0.328,0.341,0.355,0.368,0.379,0.391,0.403,0.415,0.424,0.434,0.443,0.454,0.462,0.469,0.476,0.483,0.492,0.499,0.505,0.511,0.516,0.52},
{0.,0.001,0.002,0.004,0.007,0.013,0.022,0.037,0.057,0.085,0.119,0.157,0.196,0.234,0.27,0.305,0.339,0.371,0.403,0.434,0.464,0.495,0.523,0.552,0.581,0.607,0.633,0.659,0.686,0.707,0.729,0.753,0.774,0.793,0.814,0.836,0.854,0.868,0.882,0.898,0.917,0.931,0.944,0.953,0.962,0.973,0.984,0.993,1.001,1.009,1.017},
{0.,0.002,0.006,0.013,0.026,0.047,0.078,0.118,0.166,0.223,0.286,0.352,0.419,0.483,0.544,0.603,0.661,0.717,0.772,0.826,0.876,0.928,0.976,1.023,1.069,1.11,1.151,1.189,1.231,1.264,1.298,1.334,1.361,1.388,1.416,1.444,1.467,1.486,1.507,1.53,1.555,1.576,1.594,1.603,1.612,1.623,1.638,1.652,1.664,1.674,1.682},
{0.,0.008,0.023,0.049,0.086,0.138,0.2,0.274,0.362,0.462,0.571,0.682,0.793,0.897,0.996,1.088,1.178,1.262,1.344,1.422,1.495,1.568,1.637,1.702,1.765,1.823,1.878,1.931,1.987,2.032,2.079,2.126,2.164,2.201,2.24,2.279,2.311,2.339,2.367,2.397,2.429,2.455,2.478,2.493,2.507,2.523,2.541,2.558,2.573,2.587,2.597},
{0.,0.023,0.061,0.113,0.183,0.273,0.381,0.505,0.645,0.797,0.955,1.113,1.269,1.414,1.554,1.684,1.811,1.932,2.049,2.162,2.267,2.374,2.473,2.568,2.66,2.745,2.827,2.907,2.987,3.055,3.122,3.191,3.248,3.303,3.358,3.411,3.456,3.495,3.535,3.575,3.616,3.651,3.684,3.705,3.723,3.744,3.77,3.793,3.814,3.832,3.847},
{0.,0.048,0.117,0.211,0.331,0.476,0.643,0.831,1.039,1.263,1.496,1.731,1.964,2.187,2.402,2.605,2.806,2.997,3.181,3.357,3.522,3.685,3.837,3.982,4.12,4.247,4.368,4.483,4.597,4.698,4.797,4.891,4.971,5.051,5.132,5.208,5.272,5.328,5.386,5.444,5.503,5.558,5.61,5.641,5.667,5.694,5.731,5.766,5.796,5.823,5.845},
{0.,0.085,0.202,0.355,0.551,0.79,1.07,1.388,1.738,2.107,2.487,2.868,3.239,3.593,3.934,4.253,4.564,4.861,5.144,5.415,5.672,5.92,6.153,6.376,6.591,6.788,6.977,7.155,7.327,7.488,7.643,7.791,7.918,8.043,8.168,8.281,8.38,8.473,8.566,8.655,8.741,8.82,8.894,8.946,8.992,9.037,9.089,9.139,9.18,9.215,9.244},
{0.,0.156,0.371,0.659,1.021,1.455,1.949,2.494,3.078,3.687,4.304,4.917,5.516,6.087,6.639,7.156,7.657,8.132,8.585,9.014,9.418,9.807,10.175,10.517,10.843,11.147,11.433,11.7,11.953,12.191,12.415,12.627,12.81,12.987,13.158,13.314,13.455,13.583,13.708,13.826,13.938,14.043,14.139,14.21,14.27,14.326,14.393,14.46,14.518,14.565,14.602},
{0.,0.314,0.742,1.291,1.964,2.748,3.62,4.56,5.547,6.554,7.556,8.537,9.478,10.365,11.208,11.986,12.728,13.422,14.074,14.686,15.254,15.793,16.296,16.762,17.205,17.61,17.99,18.343,18.673,18.987,19.281,19.56,19.794,20.022,20.248,20.455,20.643,20.82,20.995,21.158,21.313,21.463,21.603,21.712,21.806,21.899,22.007,22.115,22.21,22.289,22.351},
{0.,0.629,1.444,2.446,3.625,4.943,6.358,7.827,9.318,10.794,12.225,13.594,14.882,16.082,17.205,18.232,19.204,20.1,20.938,21.721,22.446,23.131,23.779,24.384,24.948,25.475,25.977,26.454,26.907,27.332,27.733,28.111,28.466,28.82,29.167,29.492,29.793,30.072,30.348,30.617,30.879,31.139,31.381,31.594,31.789,31.974,32.164,32.349,32.518,32.656,32.777}};


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;

}

float calcStecker05(float energy, float redshift){
//EBL model 3:  Stecker, Malkan, and Scully
//              Astro-ph 0510449
//		Valid for opacities  0.01 < tau < 100


double tau1, tau2, tauvalue;

int zindex=0, MAXZINDEX=9;

double zvalue[9] = {0., 0.03, 0.117, 0.2, 0.5, 1., 2., 3., 5.};

double EMIN [9] = {80., 60., 35., 25., 15., 10., 7., 5., 4.};


double coeff[9][5] = {{0., 0., 0., 0., 0.},
                     {-0.020228, 1.28458, -29.1498, 285.131, -1024.64},
                     {0.010677, -0.238895, -1.004, 54.1465, -313.486},
		     {0.0251369, -0.932664, 11.4876, -45.9286, -12.1116},
		     {-0.0221285, 1.31079, -28.2156, 264.368, -914.546},
		     {-0.175348, 8.42014, -151.421, 1209.13, -3617.51},
		     {-0.311617, 14.5034, -252.81, 1956.45, -5671.36},
		     {-0.34995, 16.0968, -277.315, 2121.16, -6077.41},
		     {-0.321182, 14.6436, -250.109, 1897.00, -5390.55}};



if (redshift < 0.){
   std::cerr<<"Invalid redshift (z<0) ..."<<std::endl;
   redshift=0.;
   } else if (redshift > 5.){
      std::cerr<<"This model is only valid for z <= 5."<<std::endl;
      redshift=5.;
      }

double x = log10(energy*1e+09);


//Find zindex
 for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;


if (energy <= EMIN[zindex])
   return 0.;

if (zindex == 0){
   tau1 = 0.;
   tau2 = coeff[1][0]*pow(x,4.)+coeff[1][1]*pow(x, 3.)+coeff[1][2]*x*x+coeff[1][3]*x+coeff[1][4];
   tauvalue =tau2*redshift/zvalue[1];
   }
   else if (zindex < MAXZINDEX-1){
   tau1 = coeff[zindex][0]*pow(x,4.)+coeff[zindex][1]*pow(x, 3.)
                        +coeff[zindex][2]*x*x+coeff[zindex][3]*x+coeff[zindex][4];
   tau2 = coeff[zindex+1][0]*pow(x,4.)+coeff[zindex+1][1]*pow(x, 3.)
                        +coeff[zindex+1][2]*x*x+coeff[zindex+1][3]*x+coeff[zindex+1][4];
   tauvalue= tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);

   } else
      tauvalue = coeff[MAXZINDEX-1][0]*pow(x,4.)+coeff[MAXZINDEX-1][1]*pow(x, 3.)
                        +coeff[MAXZINDEX-1][2]*x*x+coeff[MAXZINDEX-1][3]*x+coeff[MAXZINDEX-1][4];


return pow(10., tauvalue);

}


float calcStecker05_FE(float energy, float redshift){
//EBL model 3:  Stecker, Malkan, and Scully
//              Astro-ph 0510449
//		Valid for opacities  0.01 < tau < 100


double tau1, tau2, tauvalue;

int zindex=0, MAXZINDEX=9;

double zvalue[9] = {0., 0.03, 0.117, 0.2, 0.5, 1., 2., 3., 5.};

double EMIN [9] = {80., 60., 35., 25., 15., 10., 7., 5., 4.};

double coeff[9][5] = {{0., 0., 0., 0., 0.},
		     {-0.020753, 1.31035, -29.6157, 288.807, -1035.21},
			 {0.022352, -0.796354, 8.95845, -24.8304, -79.0409},
		     {0.0258699, -0.960562, 11.8614, -47.9214, -8.90869},
		     {0.0241367, -0.912879, 11.7893, -54.9018, 39.2521},
		     {-0.210116, 10.0006, -178.308, 1412.01, -4190.38},
		     {-0.397521, 18.3389, -316.916, 2431.84, -6991.04},
		     {-0.344304, 15.8698, -273.942, 2099.29, -6025.38},
		     {-0.28918, 13.2673, -227.968, 1739.11, -4969.32}};


if (redshift < 0.){
   std::cerr<<"Invalid redshift (z<0) ..."<<std::endl;
   redshift=0.;
   } else if (redshift > 5.){
      std::cerr<<"This model is only valid for z <= 5."<<std::endl;
      redshift=5.;
      }

double x = log10(energy*1e+09);


//Find zindex
 for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;


if (energy <= EMIN[zindex])
   return 0.;

if (zindex == 0){
   tau1 = 0.;
   tau2 = coeff[1][0]*pow(x,4.)+coeff[1][1]*pow(x, 3.)+coeff[1][2]*x*x+coeff[1][3]*x+coeff[1][4];
   tauvalue =tau2*redshift/zvalue[1];
   }
   else if (zindex < MAXZINDEX-1){
   tau1 = coeff[zindex][0]*pow(x,4.)+coeff[zindex][1]*pow(x, 3.)
                        +coeff[zindex][2]*x*x+coeff[zindex][3]*x+coeff[zindex][4];
   tau2 = coeff[zindex+1][0]*pow(x,4.)+coeff[zindex+1][1]*pow(x, 3.)
                        +coeff[zindex+1][2]*x*x+coeff[zindex+1][3]*x+coeff[zindex+1][4];
   tauvalue= tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);

   } else
      tauvalue = coeff[MAXZINDEX-1][0]*pow(x,4.)+coeff[MAXZINDEX-1][1]*pow(x, 3.)
                        +coeff[MAXZINDEX-1][2]*x*x+coeff[MAXZINDEX-1][3]*x+coeff[MAXZINDEX-1][4];


return pow(10., tauvalue);

}



float calcKneiske_extendedEnergyRange(float energy, float redshift){
/************************************************************************
This is the same EBL model defined above but with energy range extended up 
to 10^8 GeV (as provided in the paper). 
************************************************************************/


int zindex=0, eindex=-1;

float zvalue[51];

for(int i=0; i<=50; i++) zvalue[i]= 0.1*i;

float evalue[41];

for(int i=0; i<=40; i++) evalue[i]= pow(10., 0.8+0.18*i);

//Number of energy entries in the opacity table
  int MAXEINDEX = 41;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 51;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 5.){
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
       redshift=5.;
       }
if (energy >= 1e8) {
   std::cerr<<"This EBL model is only valid for E < 10^8 GeV..."<<std::endl;
   energy = 1e8;
   } else if (energy < evalue[0]) return 0.;


float tautable[41][51] = { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.002,0.003,0.005,0.007,0.009,0.011,0.014,0.018,0.02,0.021,0.026,0.028,0.03,0.033,0.036,0.038,0.041,0.044,0.044},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.003,0.007,0.012,0.018,0.026,0.031,0.041,0.047,0.056,0.064,0.075,0.082,0.093,0.101,0.11,0.122,0.131,0.137,0.142,0.156,0.161,0.17,0.176,0.182,0.188,0.197,0.198,0.2,0.209,0.217,0.221,0.219,0.224},
{0.,0.,0.,0.,0.,0.,0.,0.,0.002,0.007,0.016,0.031,0.047,0.063,0.082,0.105,0.122,0.148,0.165,0.195,0.215,0.234,0.263,0.28,0.301,0.326,0.346,0.362,0.381,0.401,0.423,0.437,0.446,0.468,0.489,0.508,0.515,0.53,0.543,0.557,0.573,0.58,0.599,0.609,0.613,0.616,0.619,0.64,0.653,0.66,0.662},
{0.,0.,0.,0.002,0.006,0.016,0.032,0.057,0.091,0.131,0.178,0.229,0.278,0.325,0.375,0.416,0.463,0.509,0.558,0.598,0.646,0.687,0.732,0.778,0.814,0.854,0.89,0.93,0.97,0.999,1.036,1.067,1.089,1.119,1.16,1.18,1.213,1.234,1.257,1.278,1.294,1.311,1.333,1.34,1.36,1.377,1.398,1.415,1.418,1.426,1.432},
{0.,0.004,0.014,0.032,0.063,0.104,0.158,0.222,0.297,0.382,0.475,0.572,0.67,0.763,0.854,0.946,1.036,1.125,1.212,1.296,1.377,1.46,1.54,1.611,1.687,1.755,1.828,1.89,1.958,2.021,2.069,2.126,2.178,2.231,2.286,2.325,2.363,2.4,2.436,2.475,2.503,2.539,2.569,2.587,2.602,2.618,2.658,2.676,2.698,2.71,2.725},
{0.,0.021,0.055,0.107,0.178,0.267,0.377,0.509,0.661,0.826,1.002,1.184,1.365,1.545,1.721,1.894,2.063,2.229,2.388,2.542,2.693,2.838,2.975,3.107,3.229,3.348,3.466,3.576,3.676,3.772,3.863,3.948,4.019,4.107,4.19,4.266,4.325,4.376,4.439,4.491,4.542,4.59,4.638,4.688,4.721,4.751,4.786,4.817,4.841,4.865,4.889},
{0.,0.053,0.132,0.244,0.394,0.584,0.815,1.086,1.388,1.711,2.048,2.39,2.729,3.053,3.371,3.666,3.958,4.242,4.513,4.775,5.02,5.26,5.492,5.709,5.916,6.113,6.308,6.483,6.652,6.811,6.969,7.119,7.25,7.379,7.505,7.618,7.72,7.822,7.921,8.01,8.098,8.179,8.251,8.307,8.363,8.414,8.468,8.517,8.564,8.599,8.624},
{0.,0.124,0.304,0.554,0.878,1.273,1.731,2.242,2.797,3.379,3.973,4.569,5.153,5.716,6.263,6.776,7.275,7.753,8.208,8.641,9.054,9.448,9.819,10.173,10.507,10.822,11.116,11.394,11.655,11.897,12.128,12.349,12.541,12.73,12.908,13.072,13.218,13.351,13.48,13.601,13.716,13.829,13.932,14.005,14.073,14.136,14.203,14.27,14.328,14.373,14.41},
{0.,0.287,0.686,1.207,1.855,2.616,3.469,4.399,5.38,6.387,7.397,8.391,9.346,10.253,11.114,11.91,12.671,13.384,14.054,14.686,15.273,15.828,16.346,16.829,17.289,17.71,18.106,18.472,18.817,19.142,19.445,19.732,19.979,20.22,20.459,20.676,20.873,21.055,21.238,21.407,21.565,21.726,21.876,21.993,22.095,22.191,22.298,22.41,22.508,22.59,22.657},
{0.,0.61,1.407,2.396,3.568,4.886,6.311,7.798,9.313,10.818,12.284,13.689,15.012,16.247,17.404,18.463,19.467,20.394,21.26,22.069,22.824,23.533,24.203,24.832,25.42,25.968,26.491,26.988,27.46,27.9,28.315,28.705,29.075,29.446,29.812,30.152,30.463,30.746,31.031,31.311,31.582,31.852,32.103,32.325,32.528,32.719,32.912,33.1,33.271,33.411,33.539},
{0.,1.106,2.469,4.072,5.875,7.808,9.812,11.834,13.836,15.781,17.653,19.431,21.106,22.66,24.122,25.461,26.731,27.924,29.039,30.095,31.097,32.055,32.969,33.842,34.678,35.476,36.251,36.999,37.715,38.402,39.069,39.714,40.331,40.941,41.539,42.108,42.637,43.134,43.62,44.088,44.538,44.972,45.374,45.738,46.074,46.384,46.686,46.974,47.237,47.463,47.659},
{0.,1.677,3.626,5.798,8.135,10.574,13.067,15.567,18.039,20.457,22.794,25.026,27.141,29.132,31.022,32.793,34.506,36.139,37.71,39.222,40.684,42.1,43.467,44.797,46.074,47.301,48.507,49.676,50.793,51.857,52.871,53.862,54.814,55.743,56.639,57.499,58.295,59.034,59.746,60.433,61.093,61.724,62.308,62.849,63.356,63.819,64.262,64.668,65.044,65.375,65.68},
{0.,2.16,4.576,7.212,10.022,12.954,15.98,19.065,22.176,25.273,28.315,31.27,34.112,36.825,39.451,41.943,44.38,46.72,48.989,51.189,53.317,55.387,57.373,59.304,61.169,62.965,64.688,66.365,67.965,69.49,70.962,72.402,73.779,75.107,76.389,77.621,78.783,79.886,80.958,81.975,82.949,83.877,84.742,85.56,86.326,87.041,87.72,88.354,88.947,89.479,89.954},
{0.,2.543,5.375,8.487,11.855,15.458,19.275,23.262,27.35,31.481,35.585,39.599,43.494,47.235,50.869,54.325,57.722,60.995,64.177,67.267,70.269,73.198,76.029,78.771,81.432,84.017,86.518,88.942,91.276,93.529,95.719,97.856,99.905,101.899,103.847,105.738,107.527,109.228,110.884,112.481,114.018,115.507,116.921,118.279,119.573,120.795,121.98,123.119,124.23,125.281,126.279},
{0.,2.967,6.367,10.207,14.469,19.112,24.092,29.343,34.783,40.338,45.898,51.38,56.728,61.883,66.922,71.747,76.51,81.139,85.672,90.105,94.451,98.726,102.899,106.995,111.007,114.953,118.829,122.686,126.504,130.297,134.099,137.956,141.807,145.679,149.638,153.663,157.724,161.838,166.035,170.336,174.695,179.081,183.519,188.018,192.487,196.929,201.357,205.748,210.043,214.168,218.122},
{0.,3.652,7.934,12.838,18.338,24.394,30.957,37.973,45.339,52.962,60.694,68.412,76.029,83.493,90.944,98.285,105.791,113.399,121.22,129.278,137.661,146.39,155.491,164.999,174.964,185.306,196.097,207.316,218.93,231.036,243.619,256.731,270.282,284.332,298.845,313.819,329.094,344.73,360.666,376.692,392.633,408.478,424.043,439.447,454.342,468.691,482.38,495.62,508.014,519.537,530.318},
{0.,4.592,10.043,16.371,23.585,31.666,40.671,50.727,61.968,74.536,88.417,103.554,119.844,137.093,155.535,174.959,195.973,218.4,242.629,268.653,296.553,326.564,358.354,391.683,426.357,462.283,499.028,536.608,574.631,613.333,652.155,691.635,730.796,769.902,809.203,847.653,885.245,922.625,958.747,993.512,1027.036,1059.378,1089.958,1118.639,1146.002,1171.957,1195.691,1217.644,1238.596,1257.107,1274.11},
{0.,5.911,13.202,22.158,33.239,47.163,65.028,88.085,117.45,153.809,196.898,245.956,299.586,356.469,416.804,479.116,545.189,613.522,685.139,758.936,835.181,913.127,992.404,1071.295,1150.584,1229.543,1306.549,1382.978,1458.044,1530.251,1601.541,1671.298,1738.967,1803.26,1865.427,1927.034,1984.922,2040.286,2093.451,2143.816,2191.08,2235.514,2278.534,2318.042,2354.268,2387.5,2418.925,2448.089,2475.51,2499.783,2521.147},
{0.,9.009,21.905,40.924,69.354,110.389,166.619,239.642,329.369,434.502,551.374,676.077,805.074,935.13,1068.127,1198.761,1331.921,1464.838,1597.544,1729.686,1861.429,1991.163,2117.745,2242.157,2363.21,2480.272,2592.943,2700.797,2805.663,2905.764,3002.492,3095.326,3181.724,3266.074,3346.168,3422.5,3494.678,3562.859,3626.573,3686.188,3743.668,3796.101,3843.547,3887.201,3929.151,3967.5,4002.189,4035.082,4062.839,4087.74,4110.684},
{0.,22.03,58.345,114.949,196.895,305.804,441.434,601.89,783.661,981.775,1189.916,1402.18,1613.917,1821.036,2026.358,2223.118,2418.201,2607.655,2792.789,2972.401,3145.981,3315.884,3478.184,3632.535,3782.126,3923.559,4057.134,4186.818,4307.988,4421.621,4529.106,4635.345,4731.974,4821.854,4909.586,4991.282,5066.075,5135.25,5202.233,5262.064,5317.433,5370.475,5420.16,5463.793,5503.402,5539.215,5571.581,5603.335,5631.482,5655.009,5674.906},
{0.,59.652,150.25,276.051,438.318,634.789,861.566,1114.402,1386.498,1671.122,1959.655,2245.073,2521.427,2786.16,3042.478,3282.021,3517.607,3739.585,3953.734,4157.298,4352.503,4536.787,4711.166,4876.395,5032.541,5180.596,5317.014,5446.549,5567.992,5681.536,5787.791,5888.478,5981.157,6067.212,6148.353,6224.281,6295.156,6358.479,6418.829,6475.29,6527.557,6573.733,6615.369,6652.971,6688.016,6719.824,6748.95,6774.913,6799.349,6821.378,6840.438},
{0.,118.494,278.927,481.595,724.009,1000.13,1303.285,1626.507,1961.35,2300.04,2633.008,2953.707,3258.802,3544.348,3817.979,4069.297,4310.603,4537.821,4752.204,4951.775,5139.664,5319.047,5484.217,5638.147,5784.602,5920.299,6045.535,6163.309,6272.021,6371.124,6465.534,6555.511,6635.531,6710.66,6781.522,6847.047,6906.352,6960.062,7010.595,7056.398,7099.781,7138.903,7174.945,7206.892,7235.127,7259.714,7283.96,7306.937,7326.664,7343.925,7359.449},
{0.,179.343,403.482,668.36,966.825,1290.621,1631.734,1983.335,2336.308,2684.14,3018.321,3334.5,3630.858,3904.786,4162.6,4397.559,4620.018,4826.435,5018.161,5197.133,5363.299,5518.294,5663.241,5796.275,5921.378,6036.505,6143.151,6241.814,6331.801,6415.154,6493.183,6567.396,6632.354,6693.691,6750.846,6803.575,6852.335,6896.752,6938.549,6974.947,7009.145,7040.074,7067.378,7092.088,7114.832,7133.886,7152.984,7171.624,7189.091,7203.26,7214.619},
{0.,223.905,486.639,780.145,1096.24,1426.756,1764.579,2103.032,2435.165,2755.93,3058.755,3342.148,3604.041,3844.103,4068.533,4268.962,4458.561,4632.823,4794.785,4942.437,5079.423,5208.149,5325.903,5434.681,5536.922,5630.042,5714.156,5793.381,5865.615,5931.717,5993.094,6052.48,6103.568,6150.938,6196.689,6239.087,6276.18,6308.829,6341.02,6369.706,6396.636,6420.688,6441.667,6460.185,6477.395,6493.233,6507.595,6521.445,6534.891,6545.57,6554.562},
{0.,240.622,508.819,796.089,1095.344,1399.919,1703.633,2001.777,2289.522,2563.372,2819.096,3056.468,3273.865,3471.321,3654.232,3817.822,3971.044,4110.555,4239.157,4357.532,4466.75,4567.834,4660.479,4745.954,4824.617,4897.16,4963.151,5024.792,5080.57,5131.459,5178.21,5223.147,5263.375,5299.957,5334.083,5365.871,5394.335,5420.211,5444.386,5465.703,5485.696,5503.714,5519.862,5533.871,5546.482,5558.239,5569.31,5579.874,5589.993,5598.566,5605.766},
{0.,230.425,477.497,734.078,995.078,1255.276,1510.209,1757.013,1992.392,2214.309,2420.055,2609.483,2782.14,2938.057,3081.465,3209.623,3328.595,3437.413,3537.528,3628.418,3711.829,3789.643,3860.361,3925.063,3985.562,4041.236,4091.024,4137.058,4179.067,4217.457,4252.827,4286.727,4316.502,4343.948,4369.604,4393.746,4415.326,4434.192,4451.867,4467.756,4482.504,4495.454,4507.267,4517.729,4527.243,4535.838,4544.152,4552.424,4559.841,4565.887,4570.89},
{0.,202.681,413.887,628.576,843.212,1053.882,1257.986,1453.662,1638.8,1812.328,1972.309,2118.488,2251.346,2370.954,2480.718,2578.164,2668.782,2751.28,2826.663,2895.299,2958.374,3016.486,3069.574,3118.214,3163.258,3203.94,3240.938,3275.439,3306.513,3335.046,3361.552,3386.303,3407.941,3428.204,3447.177,3464.437,3479.884,3493.665,3506.807,3518.197,3528.676,3537.935,3546.198,3553.652,3560.503,3566.669,3572.549,3578.374,3583.826,3587.953,3591.393},
{0.,167.672,339.009,510.649,679.928,844.437,1002.568,1153.087,1294.873,1426.964,1548.102,1658.49,1758.473,1848.227,1930.501,2003.217,2070.778,2131.868,2187.917,2238.535,2284.667,2327.384,2366.315,2401.557,2434.226,2464.24,2491.03,2515.753,2537.95,2557.996,2576.503,2594.125,2609.473,2623.305,2636.373,2648.479,2659.148,2668.31,2676.819,2684.441,2691.551,2697.816,2703.416,2708.291,2712.807,2716.778,2720.493,2724.093,2727.433,2730.024,2732.022},
{0.,132.916,267.007,399.77,529.464,654.698,774.278,887.563,993.712,1092.009,1181.721,1263.3,1336.798,1402.501,1462.408,1515.257,1563.981,1607.879,1647.603,1683.403,1716.156,1746.1,1772.958,1797.455,1819.896,1839.977,1857.87,1874.344,1888.938,1902.108,1914.228,1925.56,1935.252,1944.066,1952.271,1959.681,1966.154,1971.705,1976.94,1981.504,1985.686,1989.34,1992.627,1995.382,1997.82,1999.95,2002.034,2003.986,2005.82,2007.309,2008.568},
{0.,102.254,204.423,304.611,401.832,495.029,583.364,666.41,743.425,814.174,878.273,935.997,987.565,1033.272,1074.528,1110.612,1143.47,1172.766,1199.214,1222.729,1243.674,1262.755,1279.814,1294.948,1308.732,1321.029,1331.865,1341.678,1350.293,1357.935,1364.801,1371.326,1376.822,1381.752,1386.24,1390.301,1393.874,1396.846,1399.637,1402.039,1404.177,1406.027,1407.659,1409.011,1410.209,1411.2,1412.325,1413.441,1414.562,1415.391,1416.009},
{0.,76.87,152.822,226.56,297.319,364.232,426.76,484.591,537.404,585.214,627.839,665.663,698.986,728.122,754.067,776.4,796.516,814.238,829.861,843.546,855.744,866.633,876.154,884.605,892.252,898.906,904.695,909.881,914.314,918.308,921.931,925.277,927.992,930.471,932.79,934.87,936.627,938.046,939.402,940.546,941.605,942.458,943.183,943.817,944.361,944.873,945.458,946.067,946.674,947.127,947.48},
{0.,56.108,110.568,162.481,211.195,256.23,297.341,334.503,367.737,397.178,422.901,445.333,464.76,481.476,496.131,508.586,519.606,529.164,537.542,544.773,551.044,556.629,561.508,565.781,569.562,572.858,575.697,578.274,580.453,582.388,584.071,585.67,586.985,588.18,589.257,590.213,591.033,591.687,592.292,592.793,593.266,593.638,593.95,594.243,594.533,594.782,595.081,595.379,595.664,595.878,596.027},
{0.,38.901,75.612,109.542,140.419,168.142,192.747,214.435,233.37,249.758,263.8,275.839,286.096,294.799,302.318,308.647,314.168,318.912,322.993,326.495,329.535,332.234,334.519,336.531,338.334,339.905,341.227,342.406,343.374,344.238,345.026,345.773,346.373,346.903,347.393,347.833,348.205,348.5,348.76,348.979,349.194,349.371,349.524,349.662,349.779,349.894,350.031,350.169,350.299,350.391,350.464},
{0.,24.974,47.678,67.898,85.693,101.188,114.579,126.096,135.924,144.276,151.318,157.274,162.286,166.503,170.091,173.104,175.7,177.911,179.814,181.445,182.838,184.058,185.117,186.041,186.85,187.544,188.129,188.657,189.108,189.515,189.86,190.189,190.463,190.709,190.927,191.117,191.279,191.409,191.53,191.628,191.716,191.786,191.844,191.9,191.958,192.009,192.07,192.132,192.193,192.238,192.275},
{0.,14.631,27.4,38.373,47.738,55.683,62.394,68.053,72.808,76.795,80.122,82.909,85.236,87.179,88.824,90.199,91.378,92.379,93.237,93.966,94.593,95.148,95.615,96.018,96.379,96.69,96.956,97.197,97.393,97.566,97.717,97.863,97.984,98.09,98.185,98.268,98.339,98.397,98.45,98.491,98.53,98.562,98.591,98.618,98.643,98.666,98.694,98.721,98.746,98.765,98.781},
{0.,7.822,14.404,19.899,24.48,28.292,31.461,34.1,36.295,38.122,39.636,40.897,41.946,42.818,43.553,44.17,44.696,45.137,45.516,45.845,46.123,46.363,46.571,46.751,46.91,47.047,47.161,47.263,47.348,47.429,47.497,47.562,47.615,47.663,47.705,47.742,47.774,47.8,47.823,47.841,47.858,47.872,47.883,47.893,47.904,47.913,47.924,47.936,47.947,47.956,47.963},
{0.,3.871,7.043,9.639,11.769,13.521,14.963,16.156,17.143,17.961,18.636,19.197,19.661,20.048,20.373,20.646,20.877,21.073,21.24,21.382,21.504,21.612,21.703,21.781,21.85,21.909,21.959,22.006,22.045,22.08,22.109,22.138,22.161,22.182,22.2,22.215,22.229,22.241,22.251,22.259,22.265,22.271,22.276,22.28,22.285,22.289,22.295,22.3,22.305,22.309,22.313},
{0.,1.814,3.276,4.459,5.42,6.205,6.849,7.379,7.816,8.178,8.476,8.723,8.927,9.097,9.24,9.36,9.462,9.546,9.62,9.684,9.738,9.784,9.824,9.858,9.888,9.915,9.937,9.957,9.973,9.989,10.002,10.015,10.025,10.034,10.042,10.049,10.056,10.061,10.065,10.069,10.072,10.074,10.076,10.078,10.08,10.082,10.085,10.088,10.09,10.092,10.093},
{0.,0.823,1.479,2.007,2.433,2.781,3.065,3.299,3.491,3.65,3.781,3.889,3.979,4.054,4.117,4.169,4.214,4.252,4.284,4.311,4.335,4.356,4.373,4.388,4.401,4.413,4.422,4.431,4.439,4.445,4.451,4.456,4.461,4.465,4.468,4.471,4.474,4.476,4.478,4.479,4.48,4.481,4.482,4.483,4.484,4.485,4.486,4.487,4.488,4.488,4.489},
{0.,0.366,0.657,0.89,1.078,1.231,1.356,1.458,1.543,1.613,1.67,1.718,1.757,1.79,1.817,1.84,1.86,1.876,1.89,1.902,1.912,1.921,1.929,1.935,1.941,1.946,1.95,1.953,1.957,1.96,1.962,1.964,1.966,1.968,1.969,1.971,1.972,1.973,1.974,1.974,1.975,1.976,1.976,1.977,1.977,1.977,1.978,1.978,1.979,1.979,1.979}};


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;

}


float calcKneiske_HighUV_extendedEnergyRange(float energy, float redshift){
/************************************************************************
This is the same EBL model defined above but with energy range extended up 
to 10^8 GeV (as provided in the paper). 
************************************************************************/


int zindex=0, eindex=-1;

float zvalue[51];

for(int i=0; i<=50; i++) zvalue[i]= 0.1*i;

float evalue[44];

for(int i=0; i<=43; i++) evalue[i]= pow(10., 0.26+0.18*i);

//Number of energy entries in the opacity table
  int MAXEINDEX = 44;
//Number of redshift entries in the opacity table
  int MAXZINDEX = 51;

float tau1, tau2, tauvalue;

if(redshift < 0.){
   std::cerr<<"Invalid redshift (z < 0)..."<<std::endl;
   redshift = 0.;
   } else if (redshift > 5.){
       std::cerr<<"This EBL model is valid only for z <= 5.0"<<std::endl;
       redshift=5.;
       }
if (energy >= 1e8) {
   std::cerr<<"This EBL model is only valid for E < 10^8 GeV..."<<std::endl;
   energy = 1e8;
   } else if (energy < evalue[0]) return 0.;

float tautable[44][51] = {
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.002,0.002},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.003,0.003,0.003,0.004,0.004,0.004,0.005,0.005,0.006,0.006,0.007,0.007,0.008,0.008,0.009,0.009,0.01,0.01,0.011,0.011,0.012,0.013,0.013,0.014},
{0.,0.,0.,0.,0.,0.,0.,0.,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.003,0.004,0.005,0.006,0.007,0.009,0.01,0.011,0.013,0.014,0.016,0.017,0.019,0.021,0.023,0.025,0.027,0.029,0.031,0.033,0.036,0.038,0.04,0.043,0.045,0.047,0.05,0.052,0.054,0.056,0.058,0.061,0.063,0.064,0.066},
{0.,0.,0.,0.,0.,0.001,0.002,0.003,0.004,0.006,0.008,0.009,0.011,0.014,0.017,0.021,0.025,0.03,0.034,0.039,0.045,0.051,0.056,0.063,0.069,0.076,0.083,0.089,0.097,0.103,0.11,0.117,0.125,0.132,0.138,0.145,0.152,0.158,0.165,0.171,0.177,0.182,0.187,0.192,0.197,0.201,0.206,0.21,0.214,0.216,0.219},
{0.,0.,0.,0.001,0.002,0.004,0.006,0.01,0.014,0.021,0.031,0.043,0.056,0.071,0.086,0.102,0.117,0.133,0.15,0.167,0.184,0.201,0.218,0.234,0.251,0.267,0.283,0.298,0.314,0.328,0.341,0.355,0.368,0.379,0.391,0.403,0.415,0.424,0.434,0.443,0.454,0.462,0.469,0.476,0.483,0.492,0.499,0.505,0.511,0.516,0.52},
{0.,0.001,0.002,0.004,0.007,0.013,0.022,0.037,0.057,0.085,0.119,0.157,0.196,0.234,0.27,0.305,0.339,0.371,0.403,0.434,0.464,0.495,0.523,0.552,0.581,0.607,0.633,0.659,0.686,0.707,0.729,0.753,0.774,0.793,0.814,0.836,0.854,0.868,0.882,0.898,0.917,0.931,0.944,0.953,0.962,0.973,0.984,0.993,1.001,1.009,1.017},
{0.,0.002,0.006,0.013,0.026,0.047,0.078,0.118,0.166,0.223,0.286,0.352,0.419,0.483,0.544,0.603,0.661,0.717,0.772,0.826,0.876,0.928,0.976,1.023,1.069,1.11,1.151,1.189,1.231,1.264,1.298,1.334,1.361,1.388,1.416,1.444,1.467,1.486,1.507,1.53,1.555,1.576,1.594,1.603,1.612,1.623,1.638,1.652,1.664,1.674,1.682},
{0.,0.008,0.023,0.049,0.086,0.138,0.2,0.274,0.362,0.462,0.571,0.682,0.793,0.897,0.996,1.088,1.178,1.262,1.344,1.422,1.495,1.568,1.637,1.702,1.765,1.823,1.878,1.931,1.987,2.032,2.079,2.126,2.164,2.201,2.24,2.279,2.311,2.339,2.367,2.397,2.429,2.455,2.478,2.493,2.507,2.523,2.541,2.558,2.573,2.587,2.597},
{0.,0.023,0.061,0.113,0.183,0.273,0.381,0.505,0.645,0.797,0.955,1.113,1.269,1.414,1.554,1.684,1.811,1.932,2.049,2.162,2.267,2.374,2.473,2.568,2.66,2.745,2.827,2.907,2.987,3.055,3.122,3.191,3.248,3.303,3.358,3.411,3.456,3.495,3.535,3.575,3.616,3.651,3.684,3.705,3.723,3.744,3.77,3.793,3.814,3.832,3.847},
{0.,0.048,0.117,0.211,0.331,0.476,0.643,0.831,1.039,1.263,1.496,1.731,1.964,2.187,2.402,2.605,2.806,2.997,3.181,3.357,3.522,3.685,3.837,3.982,4.12,4.247,4.368,4.483,4.597,4.698,4.797,4.891,4.971,5.051,5.132,5.208,5.272,5.328,5.386,5.444,5.503,5.558,5.61,5.641,5.667,5.694,5.731,5.766,5.796,5.823,5.845},
{0.,0.085,0.202,0.355,0.551,0.79,1.07,1.388,1.738,2.107,2.487,2.868,3.239,3.593,3.934,4.253,4.564,4.861,5.144,5.415,5.672,5.92,6.153,6.376,6.591,6.788,6.977,7.155,7.327,7.488,7.643,7.791,7.918,8.043,8.168,8.281,8.38,8.473,8.566,8.655,8.741,8.82,8.894,8.946,8.992,9.037,9.089,9.139,9.18,9.215,9.244},
{0.,0.156,0.371,0.659,1.021,1.455,1.949,2.494,3.078,3.687,4.304,4.917,5.516,6.087,6.639,7.156,7.657,8.132,8.585,9.014,9.418,9.807,10.175,10.517,10.843,11.147,11.433,11.7,11.953,12.191,12.415,12.627,12.81,12.987,13.158,13.314,13.455,13.583,13.708,13.826,13.938,14.043,14.139,14.21,14.27,14.326,14.393,14.46,14.518,14.565,14.602},
{0.,0.314,0.742,1.291,1.964,2.748,3.62,4.56,5.547,6.554,7.556,8.537,9.478,10.365,11.208,11.986,12.728,13.422,14.074,14.686,15.254,15.793,16.296,16.762,17.205,17.61,17.99,18.343,18.673,18.987,19.281,19.56,19.794,20.022,20.248,20.455,20.643,20.82,20.995,21.158,21.313,21.463,21.603,21.712,21.806,21.899,22.007,22.115,22.21,22.289,22.351},
{0.,0.629,1.444,2.446,3.625,4.943,6.358,7.827,9.318,10.794,12.225,13.594,14.882,16.082,17.205,18.232,19.204,20.1,20.938,21.721,22.446,23.131,23.779,24.384,24.948,25.475,25.977,26.454,26.907,27.332,27.733,28.111,28.466,28.82,29.167,29.492,29.793,30.072,30.348,30.617,30.879,31.139,31.381,31.594,31.789,31.974,32.164,32.349,32.518,32.656,32.777},
{0.,1.114,2.478,4.073,5.858,7.764,9.731,11.712,13.669,15.566,17.384,19.108,20.728,22.231,23.643,24.936,26.162,27.312,28.39,29.409,30.372,31.298,32.183,33.021,33.831,34.603,35.35,36.07,36.763,37.433,38.082,38.712,39.308,39.897,40.472,41.019,41.536,42.027,42.505,42.962,43.401,43.823,44.215,44.568,44.892,45.195,45.495,45.778,46.035,46.253,46.441},
{0.,1.669,3.598,5.743,8.046,10.442,12.883,15.323,17.731,20.08,22.343,24.501,26.544,28.468,30.294,32.006,33.662,35.242,36.764,38.231,39.646,41.021,42.35,43.64,44.881,46.079,47.25,48.386,49.473,50.518,51.512,52.485,53.413,54.317,55.184,56.016,56.8,57.532,58.234,58.905,59.549,60.168,60.742,61.27,61.763,62.213,62.646,63.044,63.416,63.741,64.036},
{0.,2.136,4.521,7.115,9.873,12.744,15.697,18.701,21.726,24.732,27.681,30.541,33.292,35.923,38.47,40.89,43.259,45.534,47.745,49.889,51.961,53.981,55.92,57.803,59.624,61.381,63.062,64.698,66.257,67.76,69.207,70.616,71.956,73.257,74.511,75.704,76.846,77.933,78.988,79.983,80.934,81.848,82.702,83.503,84.249,84.938,85.598,86.22,86.807,87.325,87.788},
{0.,2.509,5.297,8.352,11.653,15.177,18.903,22.789,26.771,30.79,34.778,38.678,42.462,46.102,49.637,53.005,56.317,59.507,62.613,65.63,68.56,71.421,74.184,76.864,79.464,81.991,84.434,86.801,89.077,91.291,93.437,95.528,97.525,99.478,101.384,103.218,104.979,106.647,108.268,109.833,111.337,112.796,114.193,115.516,116.775,117.962,119.114,120.229,121.324,122.352,123.325},
{0.,2.921,6.262,10.029,14.207,18.753,23.621,28.75,34.057,39.471,44.882,50.212,55.409,60.423,65.324,70.022,74.661,79.17,83.59,87.91,92.145,96.32,100.384,104.382,108.297,112.146,115.932,119.69,123.406,127.12,130.836,134.602,138.345,142.133,146.007,149.915,153.895,157.907,162.001,166.208,170.468,174.753,179.092,183.481,187.849,192.185,196.505,200.786,204.988,209.045,212.901},
{0.,3.594,7.803,12.62,18.017,23.947,30.365,37.215,44.396,51.814,59.323,66.804,74.184,81.419,88.642,95.766,103.052,110.434,118.029,125.849,133.987,142.472,151.297,160.533,170.215,180.258,190.756,201.649,212.927,224.715,236.971,249.748,262.927,276.638,290.807,305.369,320.33,335.574,351.094,366.793,382.402,397.859,413.044,428.119,442.673,456.708,470.155,483.074,495.174,506.459,517.008},
{0.,4.522,9.883,16.1,23.173,31.079,39.871,49.668,60.592,72.768,86.177,100.749,116.399,132.973,150.711,169.421,189.695,211.341,234.768,259.933,286.965,316.065,346.878,379.241,412.896,447.814,483.551,520.094,557.154,594.81,632.683,671.247,709.468,747.585,785.933,823.557,860.285,896.809,931.969,966.084,999.003,1030.573,1060.344,1088.655,1115.352,1140.752,1164.162,1185.473,1205.809,1224.05,1240.837},
{0.,5.817,12.982,21.773,32.626,46.227,63.617,85.995,114.397,149.43,190.839,237.882,289.254,343.814,401.829,461.892,525.716,591.792,661.184,732.746,806.81,882.462,959.504,1036.319,1113.382,1190.331,1265.281,1339.8,1413.2,1483.491,1552.994,1621.504,1687.392,1750.21,1811.088,1871.193,1927.922,1982.206,2033.83,2083.211,2129.596,2172.987,2214.896,2254.042,2289.717,2322.074,2352.767,2381.021,2407.781,2431.555,2452.704},
{0.,8.862,21.536,40.189,67.983,107.946,162.51,233.157,319.749,420.996,533.425,653.333,777.362,902.585,1030.929,1157.256,1286.203,1415.058,1543.933,1672.327,1800.501,1926.675,2049.972,2171.142,2288.964,2403.087,2512.951,2618.281,2720.847,2818.268,2912.608,3003.61,3088.062,3170.439,3248.702,3323.602,3394.443,3461.089,3522.71,3581.239,3637.972,3688.986,3735.056,3778.213,3819.382,3857.135,3891.369,3923.312,3950.377,3974.814,3997.251},
{0.,21.666,57.282,112.597,192.408,298.16,429.592,584.852,760.38,951.454,1152.171,1356.795,1560.959,1760.923,1959.428,2150.012,2339.194,2523.196,2703.201,2877.832,3046.963,3212.376,3370.264,3520.875,3666.709,3804.684,3934.933,4061.468,4180.031,4290.809,4395.739,4500.048,4594.359,4682.24,4768.268,4848.114,4921.414,4989.073,5054.016,5112.905,5167.597,5219.316,5267.638,5310.852,5349.897,5384.837,5416.48,5447.072,5474.31,5497.296,5517.126},
{0.,58.399,146.872,269.436,427.219,617.897,837.646,1082.46,1345.543,1620.483,1899.194,2174.806,2441.686,2697.706,2945.96,3178.31,3406.997,3622.769,3831.204,4029.198,4219.652,4399.116,4568.833,4730.152,4882.622,5027.053,5160.069,5286.562,5405.476,5516.103,5619.862,5718.821,5809.522,5893.563,5972.97,6047.404,6116.87,6178.708,6236.928,6292.655,6344.651,6389.495,6429.92,6467.057,6501.53,6532.869,6561.466,6586.374,6609.835,6631.25,6650.012},
{0.,115.66,272.055,469.35,704.948,972.926,1266.818,1580.013,1904.007,2231.471,2553.456,2863.563,3158.689,3435.212,3700.473,3944.483,4178.917,4399.979,4608.746,4803.01,4986.442,5161.201,5321.905,5472.387,5615.503,5747.839,5869.917,5984.958,6091.521,6188.142,6280.32,6368.729,6447.044,6520.479,6589.958,6654.006,6712.001,6764.417,6813.348,6858.642,6901.79,6939.947,6975.017,7006.556,7034.364,7058.438,7082.069,7103.974,7122.729,7139.539,7155.131},
{0.,174.837,393.172,650.933,941.02,1255.37,1586.193,1927.072,2268.912,2605.59,2929.129,3235.17,3522.079,3787.594,4037.76,4266.039,4482.282,4683.203,4870.025,5044.211,5206.572,5357.594,5498.581,5628.688,5751.071,5863.358,5967.288,6063.62,6151.831,6233.136,6309.395,6382.381,6445.996,6505.877,6561.908,6613.503,6661.109,6704.401,6744.782,6780.917,6815.113,6845.216,6871.792,6896.195,6918.585,6937.281,6955.882,6973.598,6990.084,7003.756,7015.154},
{0.,218.136,473.968,759.576,1066.815,1387.787,1715.574,2043.884,2365.718,2676.377,2969.734,3244.243,3498.006,3730.843,3948.725,4143.567,4327.938,4497.627,4655.461,4799.229,4933.121,5058.532,5173.113,5279.532,5379.565,5470.358,5552.308,5629.744,5700.641,5765.073,5825.021,5883.425,5933.498,5979.751,6024.555,6065.97,6102.173,6134.032,6165.165,6193.669,6220.592,6244.033,6264.476,6282.768,6299.728,6315.255,6329.189,6342.289,6354.949,6365.266,6374.314},
{0.,234.357,495.485,775.048,1065.993,1361.873,1656.681,1946.011,2224.983,2490.369,2738.244,2968.289,3179.028,3370.623,3548.261,3707.335,3856.367,3992.25,4117.611,4232.874,4339.635,4438.133,4528.275,4611.92,4688.916,4759.651,4823.933,4884.156,4938.91,4988.539,5034.227,5078.435,5117.862,5153.576,5187.006,5218.056,5245.8,5271.03,5294.41,5315.619,5335.647,5353.213,5368.97,5382.811,5395.244,5406.771,5417.5,5427.458,5436.928,5445.195,5452.451},
{0.,224.403,464.964,714.701,968.504,1221.344,1468.873,1708.449,1936.718,2151.842,2351.331,2534.973,2702.402,2853.725,2993.032,3117.678,3233.42,3339.425,3437.016,3525.525,3607.089,3682.905,3751.718,3815.041,3874.258,3928.535,3977.042,4022.042,4063.292,4100.724,4135.289,4168.638,4197.829,4224.62,4249.743,4273.316,4294.347,4312.738,4329.833,4345.656,4360.444,4373.07,4384.6,4394.949,4404.332,4412.758,4420.807,4428.608,4435.531,4441.353,4446.411},
{0.,197.4,403.061,612.053,820.807,1025.563,1223.782,1413.762,1593.344,1761.6,1916.756,2058.502,2187.365,2303.472,2410.117,2504.909,2593.075,2673.449,2746.946,2813.784,2875.465,2932.085,2983.75,3031.361,3075.45,3115.115,3151.158,3184.887,3215.412,3243.231,3269.135,3293.488,3314.708,3334.488,3353.068,3369.918,3384.962,3398.396,3411.105,3422.463,3432.988,3442.017,3450.084,3457.468,3464.229,3470.271,3475.953,3481.427,3486.503,3490.466,3493.955},
{0.,163.307,330.155,497.254,661.914,821.827,975.418,1121.578,1259.131,1387.231,1504.734,1611.796,1708.789,1795.931,1875.877,1946.623,2012.365,2071.893,2126.542,2175.842,2220.965,2262.586,2300.475,2334.983,2366.972,2396.233,2422.328,2446.504,2468.321,2487.869,2505.955,2523.303,2538.358,2551.861,2564.661,2576.475,2586.862,2595.793,2604.018,2611.63,2618.78,2624.891,2630.363,2635.193,2639.658,2643.55,2647.123,2650.479,2653.556,2656.036,2658.087},
{0.,129.463,260.049,389.308,515.474,637.221,753.38,863.397,966.391,1061.734,1148.772,1227.91,1299.229,1363.038,1421.268,1472.701,1520.126,1562.913,1601.658,1636.536,1668.583,1697.766,1723.91,1747.906,1769.892,1789.47,1806.902,1823.018,1837.371,1850.219,1862.069,1873.226,1882.741,1891.345,1899.385,1906.616,1912.91,1918.313,1923.368,1927.944,1932.174,1935.735,1938.95,1941.68,1944.096,1946.187,1948.181,1949.981,1951.641,1953.052,1954.343},
{0.,99.601,199.102,296.651,391.232,481.839,567.655,648.315,723.055,791.7,853.914,909.939,960.01,1004.428,1044.554,1079.691,1111.692,1140.264,1166.072,1188.988,1209.499,1228.102,1244.706,1259.542,1273.064,1285.052,1295.606,1305.207,1313.689,1321.148,1327.865,1334.293,1339.697,1344.503,1348.901,1352.859,1356.319,1359.205,1361.902,1364.331,1366.518,1368.32,1369.92,1371.257,1372.446,1373.423,1374.491,1375.506,1376.503,1377.281,1377.928},
{0.,74.875,148.841,220.636,289.475,354.533,415.289,471.476,522.755,569.171,610.575,647.317,679.703,708.041,733.298,755.062,774.664,791.96,807.215,820.555,832.511,843.126,852.388,860.685,868.197,874.679,880.311,885.386,889.761,893.66,897.207,900.509,903.184,905.598,907.869,909.886,911.573,912.948,914.262,915.433,916.527,917.363,918.082,918.711,919.251,919.751,920.296,920.838,921.36,921.781,922.16},
{0.,54.646,107.679,158.223,205.619,249.418,289.38,325.504,357.795,386.404,411.413,433.224,452.123,468.395,482.672,494.819,505.565,514.897,523.08,530.13,536.282,541.724,546.472,550.671,554.388,557.595,560.353,562.88,565.037,566.925,568.571,570.151,571.449,572.611,573.662,574.587,575.368,576.003,576.59,577.109,577.601,577.968,578.282,578.575,578.864,579.106,579.382,579.643,579.881,580.079,580.244},
{0.,37.884,73.634,106.675,136.727,163.702,187.634,208.731,227.142,243.08,256.743,268.459,278.446,286.922,294.251,300.427,305.813,310.447,314.433,317.848,320.83,323.459,325.686,327.663,329.434,330.962,332.246,333.405,334.367,335.209,335.98,336.716,337.309,337.823,338.3,338.726,339.079,339.365,339.62,339.847,340.069,340.244,340.399,340.54,340.659,340.771,340.896,341.018,341.127,341.211,341.289},
{0.,24.322,46.437,66.134,83.463,98.55,111.584,122.794,132.356,140.485,147.341,153.14,158.022,162.131,165.629,168.57,171.104,173.263,175.123,176.714,178.08,179.268,180.302,181.209,182.002,182.677,183.246,183.766,184.214,184.611,184.949,185.273,185.543,185.781,185.993,186.177,186.331,186.458,186.575,186.677,186.768,186.837,186.896,186.955,187.015,187.065,187.121,187.176,187.227,187.268,187.306},
{0.,14.254,26.695,37.39,46.515,54.255,60.79,66.302,70.93,74.811,78.052,80.767,83.033,84.927,86.531,87.874,89.024,90.002,90.841,91.552,92.167,92.707,93.163,93.559,93.913,94.215,94.474,94.711,94.906,95.075,95.223,95.367,95.486,95.588,95.681,95.761,95.829,95.886,95.937,95.979,96.019,96.051,96.08,96.109,96.135,96.158,96.183,96.207,96.228,96.245,96.262},
{0.,7.624,14.039,19.398,23.863,27.578,30.665,33.235,35.373,37.152,38.627,39.856,40.878,41.727,42.444,43.047,43.56,43.991,44.361,44.682,44.955,45.189,45.391,45.569,45.724,45.857,45.969,46.068,46.153,46.232,46.299,46.363,46.415,46.461,46.502,46.539,46.569,46.595,46.616,46.635,46.652,46.666,46.677,46.688,46.699,46.708,46.719,46.729,46.739,46.747,46.754},
{0.,3.774,6.867,9.399,11.476,13.184,14.589,15.751,16.713,17.509,18.167,18.713,19.165,19.542,19.859,20.126,20.351,20.543,20.706,20.845,20.965,21.069,21.158,21.235,21.303,21.36,21.409,21.455,21.494,21.528,21.556,21.585,21.608,21.628,21.645,21.66,21.673,21.685,21.695,21.702,21.709,21.714,21.719,21.724,21.729,21.734,21.739,21.744,21.748,21.751,21.755},
{0.,1.769,3.195,4.349,5.286,6.051,6.679,7.195,7.621,7.973,8.263,8.504,8.703,8.869,9.008,9.126,9.225,9.308,9.379,9.442,9.495,9.54,9.578,9.612,9.642,9.667,9.689,9.709,9.725,9.74,9.753,9.765,9.776,9.785,9.793,9.8,9.805,9.81,9.815,9.818,9.821,9.823,9.825,9.827,9.83,9.832,9.835,9.837,9.839,9.841,9.842},
{0.,0.802,1.443,1.958,2.374,2.712,2.989,3.217,3.404,3.559,3.686,3.792,3.88,3.953,4.014,4.065,4.109,4.146,4.177,4.204,4.227,4.247,4.264,4.279,4.292,4.303,4.312,4.321,4.329,4.335,4.34,4.346,4.35,4.354,4.357,4.36,4.363,4.365,4.366,4.368,4.369,4.37,4.371,4.372,4.373,4.374,4.375,4.376,4.377,4.377,4.378},
{0.,0.357,0.641,0.868,1.051,1.2,1.322,1.422,1.505,1.572,1.628,1.675,1.713,1.745,1.772,1.794,1.813,1.829,1.843,1.855,1.865,1.873,1.881,1.887,1.892,1.897,1.901,1.905,1.908,1.911,1.913,1.916,1.918,1.919,1.921,1.922,1.923,1.924,1.925,1.926,1.926,1.927,1.927,1.928,1.928,1.929,1.929,1.929,1.93,1.93,1.93}};


  //Determine redshift index...
  for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;

  // Determine energy index
  for(int i=0; i<MAXEINDEX-1; i++)
    if(energy >= evalue[i] && energy < evalue[i+1]) eindex = i;
  if(energy >= evalue[MAXEINDEX-1]) eindex = MAXEINDEX-1;


  if (zindex < MAXZINDEX-1){
  //Find tau for redshifts above and below source by extrapolating in energy
    if(eindex < MAXEINDEX-1){
    tau1 = tautable[eindex][zindex]+(tautable[eindex+1][zindex]-tautable[eindex][zindex])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
    tau2 = tautable[eindex][zindex+1]+(tautable[eindex+1][zindex+1]-tautable[eindex][zindex+1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);    }
     else{
       tau1=tautable[MAXEINDEX-1][zindex];
       tau2=tautable[MAXEINDEX-1][zindex+1];
       }
  //  extrapolate now in redshift
  tauvalue =tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);
  } else{
      if(eindex < MAXEINDEX-1)
        tauvalue = tautable[eindex][MAXZINDEX-1]+(tautable[eindex+1][MAXZINDEX-1]-tautable[eindex][MAXZINDEX-1])*(energy-evalue[eindex])/(evalue[eindex+1]-evalue[eindex]);
       else tauvalue = tautable[MAXEINDEX-1][MAXZINDEX-1];
	}

return tauvalue;

}


} // namespace IRB
