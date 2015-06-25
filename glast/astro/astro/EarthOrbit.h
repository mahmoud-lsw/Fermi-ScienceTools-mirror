/** @file EarthOrbit.h
    @brief definition of EarthOrbit class

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/astro/astro/EarthOrbit.h,v 1.15 2007/03/15 14:27:54 burnett Exp $
*/
#ifndef astro_EarthOrbit_H
#define astro_EarthOrbit_H


#include "astro/JulianDate.h"
#include "astro/SkyDir.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
namespace astro {
    
/** @class EarthOrbit
* @brief Postition of Earth satellite
* @author G. Tosti original code 
* @author T. Burnett convert to class
    */
    class EarthOrbit   {
    public:
        
        //! specify launch date in constructor
        EarthOrbit(JulianDate launch= JulianDate(2008 ,1,1,0.0));
        
        /** 
        * set up for calculation, with orbital parameters currently wired in
        */
        void initialize();
        
        /** Orbit calculation 
        * @param jd the Julian Date for the orbit calculation
        *  @return  The orbit position (in inertial coordinates) in km
        */
        CLHEP::Hep3Vector position(JulianDate jd)const;
        
        /**
        * return the inclination of orbit
        */
        double inclination(){return s_incl;}
        
        /**
        * return the orbital phase, in terms of 'phase since ascending node was passed'
        * @param jd the Julian Date for the orbit calculation
        */
        double phase(JulianDate jd) const;
        
        /** Return the julian date
        * @param seconds number of "elapsed seconds", as exists in FluxSvc
        */ 
        ///this is for interfacing with FluxSvc, which uses "elapsed seconds" as the time parameter.
        JulianDate dateFromSeconds(double seconds) const;

        /** Return the modified julian date
         * @param seconds number of "elapsed seconds", as exists in FluxSvc
         */
        JulianDate mjdFromSeconds(double seconds) const;

        /** Return the timing correction for Shapiro delay due to the gravitational well of the sun
          * @param jd JulianDate of observation
          * @param sourceDir SkyDir giving the location of the source
          * @return Correction in seconds added to correct for the Shapiro delay
          */
        double calcShapiroDelay(JulianDate jd, const SkyDir &sourceDir) const;

        /** Return the timing correction for light travel time to solar system barycenter
          * @param jd JulianDate of observation
          * @param sourceDir SkyDir giving location of the source
          * @return Correction in seconds added to correct for light travel time
          */
        double calcTravelTime(JulianDate jd, const SkyDir &sourceDir) const;


        /** Return the correction that is added to TT to obtain TDB time
          * @param jd Julian Date in TT
          * @return Correction in seconds added to TT to obtain TDB
          */
        double tdb_minus_tt(JulianDate jd) const;

        /// set the inclination (degrees)
        static double set_inclination(double inclination);
        
    private:
        
    /**
    @brief calculate correction to phase for eccentricity 
        */
        static double Kepler(double MeanAnomaly,double Eccentricity);


        double ctatv(int long jdno, double fjdno) const;
        
        double m_M0;
        double m_dMdt;
        double m_Omega0;
        double m_dOmegadt;
        double m_w0;
        double m_dwdt;
        
        double m_a, m_alt;
        JulianDate m_launch;
        
        static double s_altitude; //<! nominal altitude (km)
        static double s_incl;     //<! orbit inclination in degrees
        static double s_e;        //<! eccentricity
        
        static double s_radius;   //<! Radius of earth
        
    };
    
}
#endif
