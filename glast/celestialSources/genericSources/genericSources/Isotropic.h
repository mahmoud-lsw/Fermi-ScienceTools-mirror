/**
 * @file Isotropic.h
 * @brief Isotropic diffuse emission.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/genericSources/genericSources/Isotropic.h,v 1.3 2012/04/24 18:21:04 jchiang Exp $
 */

#ifndef mySpectrum_Isotropic_h
#define mySpectrum_Isotropic_h

#include "flux/Spectrum.h"

/**
 * @class Isotropic
 *
 */

class Isotropic : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   Isotropic(const std::string &params);

   virtual ~Isotropic(){}

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
   virtual std::string title() const {return "Isotropic";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   double m_flux;
   double m_gamma;
   double m_emin;
   double m_emax;
   double m_ra;
   double m_dec;
   double m_cos_thetamax;

};

#endif // mySpectrum_Isotropic_h
