/** 
 * @file LaunchDirection.h
 * @brief Declare LaunchDirection class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/flux/LaunchDirection.h,v 1.12 2010/04/26 02:46:16 burnett Exp $
 */

#ifndef flux_LaunchDirection_h
#define flux_LaunchDirection_h

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/Rotation.h"


#include "astro/SkyDir.h"
#include "astro/GPS.h"

#include <algorithm>
#include <sstream>
#include <iostream>

/** @class LaunchDirection
@brief launch strategy base class

*/
class LaunchDirection  {
public:
    LaunchDirection():m_skydir(false),m_radius(0){}

    virtual ~LaunchDirection(){}

    LaunchDirection(double theta, double phi, std::string frame,double radius=0)
        :m_skydir(false)
        , m_radius(radius*M_PI/180),m_frame(frame)
    {
        CLHEP::Hep3Vector dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        m_dir = -dir;
    }
    LaunchDirection(astro::SkyDir sky, double radius=0)
        :m_skydir(true)
        , m_radius(radius*M_PI/180)
    {
        m_dir = sky.dir(); // note that this is a SkyDir, to be used as arg to LATdirection below
    }
    /** @brief choose a direction
    @param KE kinetic energy
    @param time mission time
    */
   virtual void execute(double /*KE*/, double time){
       using astro::GPS;
        if(m_skydir){
            //here, we have a SkyDir, so we need the transformation from a SkyDir to GLAST.
            m_lat_dir = GPS::instance()->LATdirection(GPS::CELESTIAL, m_dir, time); 
        }else{
            if(m_frame=="zenith"){
                //The direction is in the earth zenith system, and the rotation to GLAST is needed:
                m_lat_dir = GPS::instance()->LATdirection(GPS::GLAST, m_dir, time);
            }else{
                //otherwise, the direction is in the spacecraft system, and the rotation to GLAST is the identity:
                m_lat_dir = m_dir;
            }
        }
    }

    const CLHEP::Hep3Vector& operator()()const {return dir();}

    virtual const CLHEP::Hep3Vector& dir()const {
        static CLHEP::Hep3Vector rdir;
        rdir = m_lat_dir;
        if( m_radius>0 ) {
            // spread uniformly about a disk
            // rotate about perpendicular then about the original 
            CLHEP::Hep3Vector t(rdir);
            t.rotate( m_radius*(sqrt(CLHEP::RandFlat::shoot())),  rdir.orthogonal()),  // rotate about the orthogonal
            t.rotate( CLHEP::RandFlat::shoot( 2*M_PI ), rdir); // rotate about the original direction
            rdir = t; //replace 
        }
        return rdir;
    }

    //! solid angle: default of 1. for a point source
    virtual double solidAngle()const {
        return 1.;
    }

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << " dir" << m_dir ;
        if( m_radius>0 ) { t << " radius " << m_radius ;}
        return t.str();
    }

    /// return the cosine of the angle between the incoming direction and the earth's zenith
    virtual double zenithCosine()const{
        if(m_skydir){
            return m_dir * (astro::GPS::instance()->zenithDir()());
        }
        //if the direction is local
        return 1.0;
    }

//    virtual const astro::SkyDir & skyDirection()const { return astro::SkyDir(m_dir); }
protected:
    CLHEP::Hep3Vector m_lat_dir; ///< direction in lat frame

private:
    CLHEP::Hep3Vector m_dir;
    bool  m_skydir;
    double m_radius;
    std::string m_frame;

};

#endif // flux_LaunchDirection_h
