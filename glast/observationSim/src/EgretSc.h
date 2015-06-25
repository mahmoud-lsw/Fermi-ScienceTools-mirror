/**
 * @file EgretSc.h
 * @brief Declaration for the EGRET spacecraft object.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/src/EgretSc.h,v 1.7 2010/06/03 04:16:37 jchiang Exp $
 */

#ifndef observationSim_EgretSc_h
#define observationSim_EgretSc_h

#include <stdexcept>

#include "observationSim/Spacecraft.h"

namespace observationSim {

/**
 * @class EgretSc
 *
 * @brief Provide spacecraft attitude and orbital position information
 * for the EGRET instrument.
 *
 * This class mostly just encapsulates spacecraft data obtained from
 * the exposure history file.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/src/EgretSc.h,v 1.7 2010/06/03 04:16:37 jchiang Exp $
 */

class EgretSc : public Spacecraft {

public:

   EgretSc(astro::SkyDir &zAxis, astro::SkyDir &xAxis, 
           double earthLon, double earthLat, bool inSaa) :
      m_zAxis(zAxis), m_xAxis(xAxis), m_earthLon(earthLon),
      m_earthLat(earthLat), m_inSaa(inSaa) {}

   EgretSc(double raz, double decz, double rax, double decx,
           double earthLon, double earthLat, bool inSaa) :
      m_zAxis(astro::SkyDir(raz, decz, astro::SkyDir::EQUATORIAL)), 
      m_xAxis(astro::SkyDir(rax, decx, astro::SkyDir::EQUATORIAL)), 
      m_earthLon(earthLon), m_earthLat(earthLat), m_inSaa(inSaa) {}

   virtual ~EgretSc() {}

   virtual astro::SkyDir zAxis(double) {return m_zAxis;}
   virtual astro::SkyDir xAxis(double) {return m_xAxis;}

   virtual double EarthLon(double) {return m_earthLon;}
   virtual double EarthLat(double) {return m_earthLat;}

   virtual CLHEP::HepRotation InstrumentToCelestial(double);

   virtual bool inSaa(double) {return m_inSaa;}

   virtual void getScPosition(double time, std::vector<double> & scPosition) {
      (void)(time);
      (void)(scPosition);
      throw std::runtime_error("EgretSc::getScPosition: not implemented");
   }

   virtual void getZenith(double time, double & ra, double & dec) {
      (void)(time);
      (void)(ra);
      (void)(dec);
      throw std::runtime_error("EgretSc::getZenith: not implemented");
   }

private:

   astro::SkyDir m_zAxis;
   astro::SkyDir m_xAxis;
   double m_earthLon;
   double m_earthLat;
   bool m_inSaa;

};

} // namespace observationSim

#endif // observationSim_EgretSc_h
