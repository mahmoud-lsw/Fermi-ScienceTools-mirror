// $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/astro/src/test/test.cxx,v 1.70 2015/01/18 00:21:33 jchiang Exp $

#include <cassert>
#include <cstdlib>
#include "astro/GPS.h"
#include "astro/SolarSystem.h"
#include "astro/EarthCoordinate.h"
#include "astro/EarthOrbit.h"
#include "astro/IGRField.h"
#include "astro/JulianDate.h"
#include "astro/SkyDir.h"
#include "astro/LatProperties.h"
#include "astro/PointingTransform.h"
#include "astro/PointingHistory.h"
#include "astro/HTM.h"
#include "astro/SkyProj.h"
#include "astro/Quaternion.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "facilities/commonUtilities.h"
#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include <cstdio>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <typeinfo>

using namespace astro;
bool testSkyDir(){
    using namespace astro;
    bool ok = true;
    double l;
    for(l = 5 ; l < 355 ; l += 5.){
        for(double b = -85 ; b < 85 ; b +=5.){
            SkyDir sd4(l,b,SkyDir::GALACTIC);
            double test = fabs(l-sd4.l()) + fabs(b - sd4.b());
            //uncomment next line for verbose output:
            //std::cout << "l in = " << l << " ,b in = " << b << 
            //   " ,l out = " << sd4.l() << " ,b out = " << sd4.b() << 
            //" ,ra out = " << sd4.ra() << " ,dec out = " << sd4.dec() << std::endl;
            if(test > 1E-3){ std::cout << "error - l,b output does not match input" << std::endl; 
            ok =false;}
        }
    }
    for(double ra = 5 ; l < 355 ; l += 5.){
        for(double dec = -85 ; dec < 85 ; dec +=5.){
            SkyDir sd4(ra,dec);
            double test = fabs(ra-sd4.ra()) + fabs(dec - sd4.dec());
            //uncomment next line for verbose output:
            //std::cout << "ra in = " << ra << " ,dec in = " << dec << 
            //    " ,ra out = " << sd4.ra() << " ,dec out = " << sd4.dec() << 
            //   " ,l out = " << sd4.l() << " ,b out = " << sd4.b() << std::endl;
            if(test > 1E-3){
                std::cout << "error - ra,dec output does not match input" << std::endl; 
                ok=false;
            }
        }
    }
    //test some known locations:
    SkyDir sd5(0.,0.,SkyDir::GALACTIC);
    //std::cout <<  "X =" << sd5.dir().x() << " , Y =" << sd5.dir().y() << " , Z =" << sd5.dir().z() << std::endl;

    std::cout << "galactic center corresponds to Ra = " << sd5.ra() << " , Dec = " << sd5.dec() << std::endl;




    return ok;
}
//------------------------------------------------------------
//          test code for SkyProj
//------------------------------------------------------------
/** @brief test a single point, transforming from world to pixel and back
*/
double test_one(double lon, double lat, const astro::SkyProj& t){
    using namespace astro;
    SkyDir dir(lon, lat);
    std::pair<double, double> pixel = dir.project(t );
    //  std::cout <<  "(" << lon << ", " << lat << ") \t-> (" << pixel.first << ", " << pixel.second << ") \n";
    //std::cout <<  "\t" << lon << "\t " << lat << "\t" << pixel.first << "\t" << pixel.second << "\n";
    return SkyDir(pixel.first,pixel.second,t).difference(dir);
}

/**
 * @brief Write out a test FITS file for SkyProj
 */
std::string writeTestFile(const std::string & ctype, double * crpix, 
                          double * crval, double * cdelt) {
   std::string filename("SkyProj_test_file.fits");
   tip::IFileSvc::instance().createFile(filename);
   std::vector<long> naxes;
   naxes.push_back(1);
   naxes.push_back(1);
//   tip::IFileSvc::instance().appendImage(filename, "", naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, 
                                                            "PRIMARY");
   tip::Header & header = image->getHeader();
   header["CTYPE1"].set("RA---" + ctype);
   header["CTYPE2"].set("DEC--" + ctype);
   header["CRPIX1"].set(crpix[0]);
   header["CRPIX2"].set(crpix[1]);
   header["CRVAL1"].set(crval[0]);
   header["CRVAL2"].set(crval[1]);
   header["CDELT1"].set(cdelt[0]);
   header["CDELT2"].set(cdelt[1]);
   delete image;
   return filename;
}

#define ASSERT_EQUALS(X, Y) assert(fabs( (X - Y)/Y ) < 1e-5)

/** @brief test a group of directions
*/
bool testSkyProj(){
    using namespace astro;

    // simple test of SkyProj
    //double crpix[]={180.5,90.5},  crval[]={0,0}, cdelt[]={-1,1}; // 1-degree CAR all sky
    double crpix[]={1,1},  crval[]={-179.5,-89.5}, cdelt[]={1,1}; // 1-degree CAR all sky
    std::string ctype("CAR");
    SkyProj proj(ctype, crpix, crval, cdelt);
    if( proj.isGalactic() ) throw std::runtime_error(" wrong return from SkyProj::isGalactic");

    // write out a FITS image file using these params and create a new
    // SkyProj from the FITS header

    std::string filename(writeTestFile(ctype, crpix, crval, cdelt));
    SkyProj fileproj(filename);
    SkyProj fproj(filename,1);
    
    // create another one to verify that it is possible to have more than one

    double zcrpix[]={180.5,90.5},  zcrval[]={0,0}, zcdelt[]={-1,1}; // 1-degree CAR all sky
    SkyProj other("AIT", zcrpix, zcrval, zcdelt,0, true);
    if( !other.isGalactic() ) throw std::runtime_error(" wrong return from SkyProj::isGalactic");

    // test SkyProj::Exception
    try{
       other.pix2sph(-1000, -1000);
       // This should have thrown; if not, indicate failure.
       throw std::runtime_error("SkyProj::Exception failed");
    } catch (std::exception &) {
       // Caught as expected, so do nothing.
    }

    double delta = 20;
    for( double dec =-90+delta/2; dec<90 ; dec+= delta){
        for( double ra = 0+delta/2; ra<360; ra+=delta){
            double test = test_one(ra, dec, proj);
            if( test > 1E-8){
                std::cout << "error - SkyProj transformation failed to invert" << std::endl;
                return false;
            }

            // proj, fileproj, and fproj should operate consistently
            std::pair<double, double> coord0(proj.sph2pix(ra, dec));
            std::pair<double, double> coord1(fileproj.sph2pix(ra, dec));
            std::pair<double, double> coord2(fproj.sph2pix(ra, dec));
            
            ASSERT_EQUALS(coord0.first, coord1.first);
            ASSERT_EQUALS(coord0.second, coord1.second);

            ASSERT_EQUALS(coord0.first, coord2.first);
            ASSERT_EQUALS(coord0.second, coord2.second);
        }
    }

    // clean up
    std::remove(filename.c_str());

    return true;
}
//---------------------------------------------------------------------
void test_insideSAA() {

    //This set tests the pre-launch polygon
    //double lat_inSAA[] = {-30., -26., -20., -17., -10.,   1.,   2.,  -3.,  -8., -12., -19., -30., -28.,   0.,  -20., -25.,  -30};
    //double lon_inSAA[] = { 45.,  41.,  31.,   9., -11.1, -34., -46., -62., -79., -85., -89., -87.,  43., -34., -20., -87.9, -20};

    //double lat_notInSAA[] = {-30, -31, -26, -25, -19, -20, -16, -17, -9, -10,   2,   1,   3,   2,  -2,  -3,  -7,  -8, -11, -12, -18, -19, -31, -30};
    //double lon_notInSAA[] = { 46,  45,  42,  41,  31,  32,   9,  10,-11, -10, -34, -33, -46, -47, -62, -63, -79, -80, -85, -86, -89, -90, -87, -88};


    // test of on-orbit polygon (as of 20 Feb 2009)
    double lat_inSAA[] = {-30.0, -22.6,   2.4,   5.1,   5.2,   4.6,   0.7,  -8.6,  -9.9, -12.5, -21.7, -30.0};
    double lon_inSAA[] = { 33.9,  24.4, -18.6, -25.7, -36.0, -41.9, -58.8, -93.0, -97.4, -98.4, -92.1, -86.1}; 

    double lat_notInSAA[] = {-30.1, -22.6,   2.6,   5.2,   5.3,   4.7,    0.7,  -8.5,  -9.8, -12.5, -21.7, -30.0};
    double lon_notInSAA[] = { 33.8,  24.6, -18.6, -25.6, -36.0, -42.0,  -58.9, -93.2, -97.6, -98.6, -92.2, -86.2}; 


    // on-orbit polygon, for reference
    //static double latv[] = {-30.0, -22.6,   2.5,   5.2,   5.2,   4.6,   0.7,  -8.6,  -9.9, -12.5, -21.7, -30.0};
    //static double lonv[] = { 33.9,  24.5, -18.6, -25.7, -36.0, -42.0, -58.8, -93.1, -97.5, -98.5, -92.1, -86.1}; 

    
    bool error_found = false;
    std::cout << "\nTesting EarthCoordinate::insideSAA...";

    // Test for success
    for (unsigned int i = 0; i < sizeof(lat_inSAA)/sizeof(double); ++i) {
        astro::EarthCoordinate earthCoord; //(lat_inSAA[i], lon_inSAA[i]);
        if (!earthCoord.insideSAA(lat_inSAA[i], lon_inSAA[i]))
        {
            std::cout << std::endl << "Error: Lat/Lon " << lat_inSAA[i] << ", " << lon_inSAA[i];
            std::cout << " should return true.  Returned false.";
            error_found = true;
        }
    }

    // Test for failure
    for (unsigned int i = 0; i < sizeof(lat_notInSAA)/sizeof(double); ++i) {
        astro::EarthCoordinate earthCoord; //(lat_notInSAA[i], lon_notInSAA[i]);
        if (earthCoord.insideSAA(lat_notInSAA[i], lon_notInSAA[i]))
        {
            std::cout << std::endl << "Error: Lat/Lon " << lat_notInSAA[i] << ", " << lon_notInSAA[i];
            std::cout << " should return false.  Returned true.";
            error_found = true;
        }
    }
    if (!error_found)
        std::cout << "Done.\n" << std::endl;
    else
        throw std::runtime_error("InsideSAA test failed."); 
}

bool testJD() {
   bool passed = true;

   /// JulianDate.cxx code claims that its algorithm gives 0.00004
   /// second consistency over 2002-2020.  An accuracy better of only
   /// 0.0001 s appears to be attained with optimized builds on 32-bit
   /// machines, so we test using multiple outcomes that are within
   /// this accuracy.

   typedef std::map<double, std::vector<std::string> > Utc_vector_t;
   Utc_vector_t utc_values;

   // Test around leap second at 2009Jan01 at 00:00:00 UTC
   utc_values[252460800] = std::vector<std::string>();
   utc_values[252460800].push_back("2008-12-31T23:59:59.0000");

   utc_values[252460801] = std::vector<std::string>();
   utc_values[252460801].push_back("2009-01-01T00:00:00.0000");
   utc_values[252460801].push_back("2008-12-31T23:59:59.9999");

   utc_values[252460802] = std::vector<std::string>();
   utc_values[252460802].push_back("2009-01-01T00:00:00.0000");
   utc_values[252460802].push_back("2008-12-31T23:59:59.9999");

   utc_values[252460803] = std::vector<std::string>();
   utc_values[252460803].push_back("2009-01-01T00:00:01.0000");
   utc_values[252460803].push_back("2009-01-01T00:00:00.9999");

   utc_values[276134401.6] = std::vector<std::string>();
   utc_values[276134401.6].push_back("2009-10-01T23:59:59.6000");
   utc_values[276134401.6].push_back("2009-10-01T23:59:59.5999");

   // A random test value from previous implementation.
   utc_values[279849601.900] = std::vector<std::string>();
   utc_values[279849601.900].push_back("2009-11-13T23:59:59.9000");

   // Test around leap second introduced 2012Jul01 at 00:00:00 UTC
   // The following test value pairs were obtained from 
   // http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
   utc_values[362793601] = std::vector<std::string>();
   utc_values[362793601].push_back("2012-06-30T23:59:59.0000");
   utc_values[362793601].push_back("2012-06-30T23:59:58.9999");

   utc_values[362793602] = std::vector<std::string>();
   utc_values[362793602].push_back("2012-07-01T00:00:00.0000");
   utc_values[362793602].push_back("2012-06-30T23:59:59.9999");

   utc_values[362793603] = std::vector<std::string>();
   utc_values[362793603].push_back("2012-07-01T00:00:00.0000");
   utc_values[362793603].push_back("2012-06-30T23:59:59.9999");

   // Test around leap second introduced 2015Jul01 at 00:00:00 UTC
   // The following test value pairs were obtained from 
   // http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
   utc_values[457401602] = std::vector<std::string>();
   utc_values[457401602].push_back("2015-06-30T23:59:59.0000");
   utc_values[457401602].push_back("2015-06-30T23:59:58.9999");

   utc_values[457401603] = std::vector<std::string>();
   utc_values[457401603].push_back("2015-07-01T00:00:00.0000");
   utc_values[457401603].push_back("2015-06-30T23:59:59.9999");

   utc_values[457401604] = std::vector<std::string>();
   utc_values[457401604].push_back("2015-07-01T00:00:00.0000");
   utc_values[457401604].push_back("2015-06-30T23:59:59.9999");

   utc_values[457401605] = std::vector<std::string>();
   utc_values[457401605].push_back("2015-07-01T00:00:01.0000");
   utc_values[457401605].push_back("2015-07-01T00:00:00.9999");

   for (Utc_vector_t::const_iterator it(utc_values.begin());
        it != utc_values.end(); ++it) {
      astro::JulianDate jd 
         = astro::JulianDate(astro::JulianDate::missionStart() 
                             + it->first/86400.);
      std::string gd = jd.getGregorianDate() ;
      if (std::count(it->second.begin(), it->second.end(), gd) == 0) {
         std::cout << "Fail MET-UT test for "<< std::fixed << it->first
                   <<"\n\t got      " << gd
                   <<"\n\t expected ";
         for (size_t i(0); i < it->second.size(); i++) {
            std::cout << "\n\t\t " << it->second.at(i) << " or ";
         }
         std::cout <<  std::endl;
         passed=false;
      }
   }

   // test MET conversion
   double MET(245000000);
   JulianDate x = JulianDate::missionStart()+MET/86400.;
   std::cout << "MET = " << int(MET);
   std::cout << " is " << x.getGregorianDate() << std::endl;
   std::cout << "Mission start is: "
             << astro::JulianDate::missionStart().getGregorianDate() 
             << std::endl;

    astro::JulianDate JD2000 = astro::JulianDate(2000,1,1,12.); //2451545

    int year, year2, month, month2, day, day2;
    double utc = 12.0, utc2;
    double leap(0);
    for(year = 2008; year <= 2010; year++) { // subset

        for(month = 1; month <= 12; month++) {
            for(day = 1; day < 28; day+=5)
            {
                for(utc = 1.; utc < 24 && passed; utc+=1/3600.){
                    JD2000 = astro::JulianDate(year, month, day, utc);
                    JD2000.getGregorianDate(year2, month2, day2, utc2);

                    if(year-year2!=0 || month-month2!=0 || day-day2!=0)
                    {
                        std::cout << "Year, month, or day conversion error!" << std::endl;
                        std::cout << year << "  " << month << "  " << day << "   " << utc << std::endl;
                        std::cout << year2 << "  " << month2 << "  " << day2 << "   " << utc2 << std::endl;
                        passed = false;
                    }
                    if( 3600.*fabs(utc-utc2) > 50e-6 ) // check to 50 usec
                    {
                        std::cout << "Time error!" << 3600*fabs(utc-utc2) << " sec for "
                            << JD2000.getGregorianDate()  << std::endl;
                        passed = false;
                    }
                    
                }
            }
        }
    }
    // quick test of leap second(s) - 2 after 2008
    JD2000 = astro::JulianDate(astro::JulianDate::missionStart()+252460802.00003/86400. );
    std::string td( JD2000.getGregorianDate() ), expect("2009-01-01T00:00:00.0000");
    if( td !=expect ){
        std::cout << "Fail leap second test: expected "<< expect <<", got " << td <<  std::endl;
        passed=false;
    }
 
    if(passed)
    {
        std::cout << "JD Conversions passed!" << std::endl;
    }
    return passed;
}
bool testHTM()
{
    using astro::HTM;
    size_t maxlevel = 6;
    astro::HTM h(maxlevel);
    //      h.dump(std::cout);
    for( size_t level = 0; level <= maxlevel; ++level) {
        double area = 0;
        for( HTM::const_iterator it = h.begin(level); it!=h.end(level); ++it){
            const HTM::Node & n = *it;
            area += n.area();
        }
        double check =  fabs(area/(4*M_PI) -1);
        if( check >1e-8) { 
            std::cout << "HTM area did not add up"<< std::endl; return false; 
        }
    }
    std::cout << "HTM check OK" << std::endl; 
    return true;
}

void checkdir(const astro::SkyDir& dir1, const astro::SkyDir& dir2) {
        assert(dir1.difference(dir2) < 1e-5); //readback should correspond to contents of FITS
}
void checkdir(double ra1, double dec1, double ra2, double dec2) {
    checkdir( astro::SkyDir(ra1, dec1), astro::SkyDir(ra2, dec2) );
}

bool test_GPS_readFitsData() {
    using astro::SkyDir;
    astro::GPS * gps = astro::GPS::instance();

    //std::string filename("../src/test/test_FT2.fits");
    std::string filename=facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("astro"), "test_FT2.fits");
    gps->setPointingHistoryFile(filename);

    double time(gps->history().startTime());
    gps->time(time);

    double lat(gps->lat()), lon(gps->lon());

    // these numbers extracted from a previous run, or the data itself. 
    ASSERT_EQUALS(gps->lat(), 28.69208); 
    ASSERT_EQUALS(gps->lon(), -91.25456);
    checkdir(gps->xAxisDir(), SkyDir( 99.46017, 0));
    checkdir(gps->zAxisDir(), SkyDir(9.460165, 63.5));
    checkdir(gps->zenithDir(), SkyDir(9.460165, 28.5));


    // this should work--testing capability to merge different FT2 files
    PointingHistory history(filename);
    history.readFitsData(filename);

    return true;
}

bool test_IGRF() {
   std::ofstream output("test_IGRF_output.txt");
   EarthOrbit earthOrbit;
   for (int year(1985); year < 2020; year++) {
      for (int month(1); month < 13; month++) {
         JulianDate jd(year, month, 15, 0);
         double met(jd.seconds() - JulianDate::missionStart().seconds());
         EarthCoordinate earthCoord(earthOrbit.position(jd), met);
         output << month << "  "
                << year << "  "
                << met << "  "
                << earthCoord.L() << "  " 
                << earthCoord.B() << "  " 
                << earthCoord.lambda() << "  ";
// The following returns zeros for some reason.
//                 << earthCoord.magnetic_field().x() << "  "
//                 << earthCoord.magnetic_field().y() << "  "
//                 << earthCoord.magnetic_field().z() << "  ";
         double longitude(earthCoord.longitude());
         double latitude(earthCoord.latitude());
         double altitude(earthCoord.altitude());
// This line is from EarthCoordinate::computeMagCoords():
         double approx_year(2001 + met/3.15e7); // year does not need to be very precise
// The "year" argument is used for interpolation in igrf_sub.cxx, so
// it does need to be precise!  Here we use this approximation to
// ensure that R() and verticalRigidityCutoff() are computed at the
// same point in time as in EarthCoordinate::computeMagCoords(). That
// code has to be fixed eventually.
         IGRField::Model().compute(latitude, longitude, altitude, approx_year);

         output 
            // << "       " 
            // << IGRField::Model().L() << "  "
            // << IGRField::Model().B() << "  "
            // << IGRField::Model().lambda()*180./M_PI << "  "
            << IGRField::Model().bEast() << "  "
            << IGRField::Model().bNorth() << "  "
            << -IGRField::Model().bDown() << "  "
            << IGRField::Model().R() << "  "
            << IGRField::Model().verticalRigidityCutoff() << std::endl;
      }
   }
   output.close();

// Test for year request outside of the valid range for IGRF-11.
   JulianDate jd(2025, 1, 1, 0);
   double met(jd.seconds() - JulianDate::missionStart().seconds());
   EarthCoordinate earthCoord(earthOrbit.position(jd), met);
   try {
      IGRField::Model().compute(earthCoord.longitude(), earthCoord.latitude(),
                                earthCoord.altitude(), 2025.);
      throw std::runtime_error("Expected exception not thrown.");
   } catch(std::runtime_error & eObj) {
      std::string message(eObj.what());
      if (message.find("Requested year, 2025, is outside "
                       "the valid range of 1900-2020") == std::string::npos) {
         throw;
      }
   }
   return true;
}

bool test_LatProperties() {
   astro::GPS * gps(astro::GPS::instance());
   std::string data_path(facilities::commonUtilities::getDataPath("astro"));
   std::string filename(facilities::commonUtilities::joinPath(data_path, 
                                                              "test_FT2.fits"));
   gps->setPointingHistoryFile(filename);
   for (size_t i(1); i < 11; i++) {
      PointingHistory history(gps->history());
      double met(i*800);
      double tstart(static_cast<int>(met/30.)*30.);
      double tstop(tstart + 30.);
      assert(0 == history(met).latProperties().lat_mode());
      assert(1 == history(met).latProperties().lat_config());
      assert(1 == history(met).latProperties().data_qual());
      if (history(met).latProperties().in_saa()) {
         assert(0 == history(met).latProperties().livetime());
      } else {
         ASSERT_EQUALS(27., history(met).latProperties().livetime());
      }
      ASSERT_EQUALS(tstart, history(met).latProperties().interval_start());
      ASSERT_EQUALS(tstop, history(met).latProperties().interval_stop());
      assert(0 == history(met).latProperties().rock_angle());
   }
   return true;
}

int main(){

    using namespace astro;
    using namespace std;
    int rc = 0;

    facilities::commonUtilities::setupEnvironment();

    try {
        if( Quaternion::test()!=0) {
            std::cerr << "Failed quaternion test" << std::endl;
            rc=1;
        }else{        std::cerr << "Quaternion tests ok" << std::endl;  }

        if( GPS::test()!=0 ) {
            std::cerr << "Failed GPS test" << std::endl;
            rc=1;
        }else{         std::cerr << "GPS tests OK" << std::endl; }


        // One needs the test data to run this.
        if (!test_GPS_readFitsData()) return 1;

        JulianDate start = JulianDate::missionStart(); 
        { // code associated with JIRA
            std::cout << "JIRA OBS-14\n";
            double secsperday(86400.); 
            astro::JulianDate mission_start(2001, 1, 1, 0); 
 
            //double t0(252460741.00005 ); //version to check insertion of leap
            //double t0(255139200.00005); // slight change to round off seconds
            double t0(252460801.00005); // display 2008 leap second

            for (double i=t0-5; i < t0+5; i++) { 
                astro::JulianDate now(mission_start + i/secsperday); 
                std::cout << static_cast<int>(i+0.5) << " " << now.getGregorianDate() << std::endl; 
            } 
        }

        double test=0;

        test += fabs(start -2451910.5);

        double ra=30,dec=50;
        SkyDir sd(ra, dec);

        double l=10.,b=-10.;
        SkyDir sd3(l,b,SkyDir::GALACTIC);

        test += fabs(l-sd3.l()) + fabs(b - sd3.b());
        test += fabs(ra-sd.ra()) +fabs(dec-sd.dec());

        EarthOrbit abcd;
        double juliandate = abcd.dateFromSeconds(0.0);


        std::cout << "EarthCoordinate:\n";
        EarthCoordinate xyza(abcd.position(juliandate),juliandate);

        std::cout << std::setprecision(6) << "\tlatitude at t0 = " << xyza.latitude()
            << " , longitude at t0 = " << xyza.longitude() << std::endl;

#if 0 // this test no longer makes sense.
        double lat=25, lon=45;
        std::cout << "EarthCoordinate("<<lat<<","<<lon<<")\n";
        EarthCoordinate ec(lat, lon);

        test += fabs(lat-ec.latitude()) + fabs(lon-ec.longitude());
//        double L = ec.L(), B=ec.B();
        std::cout << "\tL=\t"<< ec.L() 
            << "\n\tB=\t"<< ec.B() 
            << "\n\tgeolat=\t"<< ec.geolat() 
            << "\n\tgeolon=\t"<< ec.geolon() 
            << std::endl;
#endif
        // test the SkyDir difference function
        SkyDir sd2(ra+1,dec);
        double 
            diff = sd2.difference(sd);
        test+= fabs( diff/cos(dec*M_PI/180) - M_PI/180. );

        //now test the galactic transformation function:
        SkyDir zenith(20,0,astro::SkyDir::GALACTIC);
        SkyDir xdir(-70,0,astro::SkyDir::GALACTIC);
        PointingTransform trans(zenith,xdir);
        CLHEP::Hep3Vector vertical(0,0,1);
        //double templ=trans.gDir(vertical).l();	
        //double tempb=trans.gDir(vertical).b();
        test += trans.gDir(vertical).l()-20.0;
        test += trans.gDir(vertical).b();

        // make sure can set ephemeris
        JulianDate jdbary(2454101.5001062094);
        JulianDate mission(2451910.5);
        SolarSystem ss(SolarSystem::EARTH);
        CLHEP::Hep3Vector bary = ss.getBarycenter(jdbary);
        CLHEP::Hep3Vector solar = ss.getSolarVector(mission);

        //expected direction 1/1/2001 : 00:00:00
        CLHEP::Hep3Vector svec = SkyDir(18.771666*15,-23.013055)();
        test += (1-svec.dot(solar.unit()))*1e3;
        solar = ss.getSolarVector(mission+1);

        //expected direction 1/2/2001 : 00:00:00
        svec = SkyDir(18.845277*15,-22.927777)();
        test += (1-svec.dot(solar.unit()))*1e3;
        solar = ss.getSolarVector(mission+365);

        //expected direction 1/1/2002 : 00:00:00
        svec = SkyDir(18.753888*15,-23.0325)();
        test += (1-svec.dot(solar.unit()))*1e3;

        std::cout << "Barycenter coords for JD: "
            << std::setprecision(8)<< jdbary << ": "  << bary << std::endl;
        // expected, from comparison with glbary.
        CLHEP::Hep3Vector expect (84.7056619, -445.132338, -192.922157);
        test += (expect-bary).mag() *1e3 ;


        /// @todo: some test of the barycenter calculation.
        /*
        // test projection (use default)

        std::pair<double,double> proj= sd.project();
        SkyDir sd4(proj.first, proj.second, astro::SkyDir::PROJECTION);
        //THB double sd4_ra= sd4.ra(), sd4_dec=sd4.dec();
        test += sd4.difference(sd);

        */

        if( fabs(test) < 1e-3 ) {
            cout << "tests ok " << endl;

        } else {
            cout << "failed a test" << endl;
            cout << "Mission start" << start << endl;  
            cout << "SkyDir("<<ra<<","<<dec<<") " << sd.ra() << ", " << sd.dec()   << endl;
            cout << "SkyDir3("<<l<<","<<b<<") " << sd3.l() << ", " << sd3.b()   << endl;

            // run the sun and moon
            rc= 1;
        }

        // other tests
        test_insideSAA();

        if( !testSkyDir() )rc=1;

        if( !testHTM() ) rc= 1;

        if(! testSkyProj() ) rc= 1;

        if( ! testJD()) rc=1;

        if (!test_IGRF()) {
           rc = 1;
        }
        if (!test_LatProperties()) {
           rc = 1;
        }
    }catch( const std::exception& e){
        std::cerr << "Failed test because caught " <<typeid(e).name()<<" \""  
            << e.what() << "\"" << std::endl;
        rc= 1;
    }
    return rc;
}

