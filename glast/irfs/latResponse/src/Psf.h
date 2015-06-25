/**
 * @file Psf.h
 * @brief Psf class declaration.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/Psf.h,v 1.1.1.5 2011/03/20 19:24:50 elwinter Exp $
 */

#ifndef latResponse_Psf_h
#define latResponse_Psf_h

#include <cmath>

#include <map>
#include <string>
#include <vector>

#include "irfInterface/IPsf.h"
#include "irfInterface/AcceptanceCone.h"

#include "latResponse/ParTables.h"

#include "PsfBase.h"

namespace latResponse {

class PsfIntegralCache;

/**
 * @class Psf
 *
 * @brief A LAT point-spread function class for post-handoff review IRFs.
 *
 */

class Psf : public PsfBase {

public:

   Psf(const std::string & fitsfile, bool isFront=true,
       const std::string & extname="RPSF", size_t nrow=0);

   Psf(const Psf & rhs);

   virtual ~Psf();

   /// A member function returning the point-spread function value.
   /// @param appDir Apparent (reconstructed) photon direction.
   /// @param energy True photon energy in MeV.
   /// @param srcDir True photon direction.
   /// @param scZAxis Spacecraft z-axis.
   /// @param scXAxis Spacecraft x-axis.
   /// @param time Photon arrival time (MET s)
   virtual double value(const astro::SkyDir & appDir, 
                        double energy, 
                        const astro::SkyDir & srcDir, 
                        const astro::SkyDir & scZAxis,
                        const astro::SkyDir & scXAxis, 
                        double time=0) const;

   /// Return the psf as a function of instrument coordinates.
   /// @param separation Angle between apparent and true photon directions
   ///        (degrees).
   /// @param energy True photon energy (MeV).
   /// @param theta True photon inclination angle (degrees).
   /// @param phi True photon azimuthal angle measured wrt the instrument
   ///            X-axis (degrees).
   /// @param time Photon arrival time (MET s)
   virtual double value(double separation, double energy, double theta,
                        double phi, double time=0) const;


   typedef std::vector<irfInterface::AcceptanceCone *> AcceptanceConeVector_t;

   /// Angular integral of the PSF over the intersection of acceptance
   /// cones.
   virtual double 
   angularIntegral(double energy,
                   const astro::SkyDir & srcDir,
                   const astro::SkyDir & scZAxis,
                   const astro::SkyDir & scXAxis,
                   const AcceptanceConeVector_t & acceptanceCones,
                   double time=0);

   virtual double 
   angularIntegral(double energy,
                   const astro::SkyDir & srcDir,
                   double theta, 
                   double phi, 
                   const AcceptanceConeVector_t & acceptanceCones, 
                   double time=0);

   virtual double angularIntegral(double energy, double theta, double phi,
                                  double radius, double time=0) const;

   virtual Psf * clone() {
      return new Psf(*this);
   }

   /// Functions from handoff_response.
   static double old_base_function(double u, double sigma, double gamma);

   static double old_base_integral(double u, double sigma, double gamma);

   static double old_integral(double sep, double * pars);

   static double old_function(double sep, double * pars);

protected:

   /// Disable this.
   Psf & operator=(const Psf &) {
      return *this;
   }

   ParTables m_parTables;

   // PSF scaling parameters
   double m_par0;
   double m_par1;
   double m_index;

   // store all of the PSF parameters
   std::vector<double> m_psf_pars;

   mutable double m_loge_last;
   mutable double m_costh_last;

   PsfIntegralCache * m_integralCache;

   mutable double m_pars[6];

   /// Hard-wired cut-off value of scaled deviation squared used by
   /// handoff_response. This should be a parameter passed in the IRF
   /// FITS header.
   static double s_ub;

   void readScaling(const std::string & fitsfile, bool isFront,
                    const std::string & extname="PSF_SCALING_PARAMS");

private:

   double * pars(double energy, double costh) const;

   /**
    * @class PsfIntegrand
    * @brief Functor used for integrating the PSF to get the proper
    * normalization.
    */
   class PsfIntegrand {
   public:
      PsfIntegrand(double * pars) : m_pars(pars) {}

      /// @param sep angle between true direction and measured (radians)
      double operator()(double sep) const {
         return old_function(sep, m_pars)*std::sin(sep);
      }
   private:
      double * m_pars;
   };

};

} // namespace latResponse

#endif // latResponse_Psf_h
