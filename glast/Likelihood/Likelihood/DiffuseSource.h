/** 
 * @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/DiffuseSource.h,v 1.53 2014/04/22 04:56:38 jchiang Exp $
 */

#ifndef Likelihood_DiffuseSource_h
#define Likelihood_DiffuseSource_h

#include <cmath>

#include <stdexcept>

#include "optimizers/dArg.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Source.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/TrapQuad.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

   class Event;
   class MapBase;
   class Observation;

/** 
 * @class DiffuseSource
 *
 * @brief Non-point-like gamma-ray source.  The extragalactic and
 * Galactic diffuse emission components are the obvious examples, but
 * this class also includes discrete diffuse sources such as the LMC
 * or supernova remnants.
 *
 * This representation assumes that a single spectral model describes
 * the emission over the entire angular extent of the source.
 * Therefore, in order to have spectral variations across a source, it
 * must comprise sub-components that can be represented using this
 * class.
 *
 * Note that the member function double spatialDist(astro::SkyDir &)
 * returns the spatial distribution of the emission as a function of
 * direction.
 *
 */

class DiffuseSource : public Source {

public:

   /// @param spatialDist A Function object describing the spatial
   ///        distribution of the emission.
   /// @param requireExposure If true, then an ExposureMap must have
   ///        been read in so that the spatially integrated for the
   ///        spectral calculation; if false, then the map is not
   ///        integrated.
   DiffuseSource(optimizers::Function * spatialDist,
                 const Observation & observation,
                 bool requireExposure=true,
                 bool mapBasedIntegral=false);
   
   DiffuseSource(const DiffuseSource &rhs);

   virtual ~DiffuseSource() {
      delete m_spatialDist;
   }

   /// Returns photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event & evt,
			      CachedResponse * cResp = 0) const;

   /// Returns the derivative wrt to the named Parameter
   virtual double fluxDensityDeriv(const Event & evt, 
                                   const std::string & paramName,
				   CachedResponse * cResp = 0) const;

   virtual double fluxDensity(double inclination, double phi, double energy,
                              const astro::SkyDir &appDir, int evtType,
                              double time, CachedResponse* cResp = 0) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      (void)(time);
      (void)(cResp);
      return 0;
   }

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy,
                                   const astro::SkyDir &appDir,
                                   int evtType,
                                   double time,
                                   const std::string & paramName,
				   CachedResponse* cResp = 0) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      (void)(time);
      (void)(paramName);
      (void)(cResp);
      return 0;
   }

   /// Return the spatial distribution of the gamma-ray emission
   double spatialDist(const astro::SkyDir & dir) const {
      SkyDirArg SDarg(dir);
      return (*m_spatialDist)(SDarg);
   }

#ifndef SWIG
   double spatialDist(SkyDirArg dir) const {
      return (*m_spatialDist)(dir);
   }
#endif

   const optimizers::Function * spatialDist() const {
      return m_spatialDist;
   }

   virtual Source *clone() const {
      return new DiffuseSource(*this);
   }

   virtual const std::vector<double> & exposure() const {
      return m_exposure;
   }

   /// @return Photon flux integrated over the ROI energy bounds. 
   /// Units are #/cm^2/s
   virtual double flux() const;

   /// @return Derivative of integrated photon flux wrt the named parameter
   double fluxDeriv(const std::string & parName) const;

   /// @return Photon flux integrated over the given energy range.
   /// Units are #/cm^2/s
   double flux(double emin, double emax, size_t npts=100) const;

   /// @return Derivative of integrated photon flux wrt the named parameter
   /// over the given energy range.
   double fluxDeriv(const std::string & parName, 
                    double emin, double emax, size_t npts=100) const;

   /// @return Energy flux integrated over the ROI energy bounds. 
   /// Units are MeV/cm^2/s
   virtual double energyFlux() const;

   virtual void computeExposure(const std::vector<double> & energies,
                                bool verbose=false);

   /// @return Derivative of integrated energy flux wrt the named parameter
   double energyFluxDeriv(const std::string & parName) const;

   /// @return Energy flux integrated over the given energy range.
   /// Units are MeV/cm^2/s
   double energyFlux(double emin, double emax, size_t npts=100) const;

   /// @return Derivative of integrated energy flux wrt the named parameter
   /// over the given energy range.
   double energyFluxDeriv(const std::string & parName, 
                          double emin, double emax, size_t npts=100) const;

   const MapBase * mapBaseObject() const;

   MapBase * mapBaseObject();

   double angularIntegral(double energy) const;

   double diffuseResponse(const Event & evt) const;

   bool mapBasedIntegral() const;

private:

   /// spatial model
   optimizers::Function * m_spatialDist;

   bool m_mapBasedIntegral;

   template<typename Functor>
   double computeEnergyIntegral(const Functor & func, 
                                double emin, double emax, size_t npts) const {
      std::vector<double> energies;
      energies.reserve(npts);
      double estep(std::log(emax/emin)/float(npts-1));
      for (size_t k(0); k < npts; k++) {
         energies.push_back(emin*std::exp(estep*k));
      }
      return computeEnergyIntegral(func, energies);
   }

   template<typename Functor>
   double computeEnergyIntegral(const Functor & func, 
                                const std::vector<double> & energies) const {
      std::vector<double> integrand;
      integrand.reserve(energies.size());
      for (size_t k(0); k < energies.size(); k++) {
         optimizers::dArg arg(energies.at(k));
         integrand.push_back(func(arg)*angularIntegral(energies.at(k)));
      }
      bool useLog;
      TrapQuad fluxIntegral(energies, integrand, useLog=true);
      return fluxIntegral.integral();
   }

   // This function fills the Source::m_exposure array which contains
   // the integral over the source region of the spatial distribution
   // of the source times the unbinned exposure map.
   void integrateSpatialDist();

   /// @return Approximate energy-dependent outer angle (in radians) 
   ///         for diffuse response calculation.
   double psfRange(double energy) const;
};

} //namespace Likelihood

#endif // Likelihood_DiffuseSource_h