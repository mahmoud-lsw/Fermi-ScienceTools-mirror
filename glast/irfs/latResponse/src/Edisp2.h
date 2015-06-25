/**
 * @file Edisp2.h
 * @brief Class definition for Riccardo's second generation energy
 * dispersion representation
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/Edisp2.h,v 1.1.1.5.2.4 2015/04/26 06:53:49 jasercio Exp $
 */

#ifndef latResponse_Edisp2_h
#define latResponse_Edisp2_h

#include <map>
#include <vector>

#include "irfInterface/IEdisp.h"
#include "latResponse/ParTables.h"
#include "EdispInterpolator.h"

namespace latResponse {

/**
 * @class Edisp2
 * @brief Edisp2 class for Riccardo's second generation energy
 * dispersion representation
 */

class Edisp2 : public irfInterface::IEdisp {

public:

   Edisp2(const std::string & fitsfile,
          const std::string & extname="ENERGY DISPERSION",
          size_t nrow=0);

   virtual ~Edisp2();

   /// A member function returning the energy dispersion function.
   /// @param appEnergy measured photon energy in MeV.
   /// @param energy True photon energy in MeV.
   /// @param srcDir True photon direction.
   /// @param scZAxis Spacecraft z-axis.
   /// @param scXAxis Spacecraft x-axis.
   /// @param time   MET
   virtual double value(double appEnergy, 
                        double energy,
                        const astro::SkyDir & srcDir, 
                        const astro::SkyDir & scZAxis,
                        const astro::SkyDir & scXAxis,
                        double time=0) const;

   virtual double value(double appEnergy, double energy,
                        double theta, double phi,
                        double time=0) const;

   virtual irfInterface::IEdisp * clone() {
      return new Edisp2(*this);
   }

   double scaleFactor(double energy, double costheta) const;

   double evaluate(double emeas, double energy,
                   double theta, double phi, double time, 
                   double * pars) const;

   void renormalize(double logE, double costh, double * params) const;

private:

   ParTables m_parTables;

   mutable double m_loge_last;
   mutable double m_costh_last;

   mutable double m_pars[10];

   bool m_renormalized;

   std::string m_fitsfile;
   std::string m_extname;
   size_t m_nrow;

   mutable EdispInterpolator * m_interpolator;

   double * pars(double energy, double costh) const;

   std::vector<double> m_scalePars;
   double m_p1;
   double m_p2;
   double m_t0;

   double old_function(double xx, double * pars) const;

   void readScaling(const std::string & fitsfile,
                    const std::string & extname="EDISP_SCALING_PARAMS");

   void renormalize();

   class EdispIntegrand {
   public:
      EdispIntegrand(double * pars, double etrue, double scaleFactor,
                     const Edisp2 & self) 
         : m_pars(pars), m_etrue(etrue), m_scaleFactor(scaleFactor),
           m_self(self) {}

      double operator()(double emeas) const {
         double xx((emeas - m_etrue)/m_etrue/m_scaleFactor);
         return m_self.old_function(xx, m_pars)/m_etrue/m_scaleFactor;
      }

   private:
      double * m_pars;
      double m_etrue;
      double m_scaleFactor;
      const Edisp2 & m_self;
   };

};

} // namespace latResponse

#endif // latResponse_Edisp2_h
