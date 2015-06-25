/**
 * @file Spacecraft.h
 * @brief Declaration of Spacecraft base class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/observationSim/Spacecraft.h,v 1.9 2010/06/03 04:16:37 jchiang Exp $
 */

#ifndef observationSim_Spacecraft_h
#define observationSim_Spacecraft_h

#include <vector>
#include "astro/SkyDir.h"

namespace observationSim {

/**
 * @class Spacecraft
 *
 * @brief An abstract base class that defines the interface for
 * objects that provide information on spacecraft position and
 * attitude.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/observationSim/Spacecraft.h,v 1.9 2010/06/03 04:16:37 jchiang Exp $
 */

class Spacecraft {

public:

   Spacecraft() : m_livetimeFrac(1) {}

   virtual ~Spacecraft() {}

   /// Spacecraft z-axis in J2000 coordinates.
   virtual astro::SkyDir zAxis(double time) = 0;

   /// Spacecraft x-axis in J2000 coordinates.
   virtual astro::SkyDir xAxis(double time) = 0;

   /// Earth longitude in degrees.
   virtual double EarthLon(double time) = 0;

   /// Earth latitude in degrees.
   virtual double EarthLat(double time) = 0;

   /// Rotation matrix from instrument to J2000 coordinates
   virtual CLHEP::HepRotation InstrumentToCelestial(double time) = 0;

   /// true if in SAA
   virtual bool inSaa(double time) = 0;

   /// Spacecraft position in geocentric coordinates (km)
   virtual void getScPosition(double time,
                              std::vector<double> & scPosition) = 0;

   virtual void getZenith(double time, double & ra, double & dec) = 0;

   virtual double livetimeFrac(double time) const {
      (void)(time);
      return m_livetimeFrac;
   }

   virtual void setLivetimeFrac(double livetimeFrac) {
      m_livetimeFrac = livetimeFrac;
   }

private:

   double m_livetimeFrac;

};

} // namespace observationSim

#endif // observationSim_Spacecraft_h
