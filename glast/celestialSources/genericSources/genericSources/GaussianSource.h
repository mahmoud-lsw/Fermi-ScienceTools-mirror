/**
 * @file GaussianSource.h
 * @brief Extended source modeled using a 2D Gaussian.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/genericSources/genericSources/GaussianSource.h,v 1.3 2012/04/24 18:21:04 jchiang Exp $
 */

#ifndef fluxSources_GaussianSource_h
#define fluxSources_GaussianSource_h

#include "flux/Spectrum.h"

/**
 * @class GaussianSource
 * @brief 2D Gaussian extended source
 *
 */

class GaussianSource : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   GaussianSource(const std::string &params);

   virtual ~GaussianSource(){}

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const;

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Total flux (photons/m^2).
   /// @param time Simulation time in seconds.
   virtual double flux(double time) const;

   /// @return "Effective" solid angle (sr).
   virtual double solidAngle() const;

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "GaussianSource";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   double m_flux;
   double m_gamma;
   double m_major;
   double m_minor;
   double m_posAngle;
   double m_emin;
   double m_emax;

   CLHEP::HepRotation m_rot;
   CLHEP::Hep3Vector m_rotatedSrcVec;

};

#endif // fluxSources_GaussianSource_h
