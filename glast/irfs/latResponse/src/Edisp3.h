/**
 * @file Edisp3.h
 * @brief Class definition for Thibaut Desgardin's third generation energy
 * dispersion representation
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/Edisp3.h,v 1.1.2.4 2015/04/26 06:53:49 jasercio Exp $
 */

#ifndef latResponse_Edisp3_h
#define latResponse_Edisp3_h

#include <map>
#include <vector>

#include "irfInterface/IEdisp.h"
#include "latResponse/ParTables.h"
#include "EdispInterpolator.h"

namespace irfUtil {
   class IrfHdus;
}

namespace latResponse {

/**
 * @class Edisp3
 * @brief Edisp3 Class for Thibaut Desgardin's third generation energy
 * dispersion representation
 */

class Edisp3 : public irfInterface::IEdisp {

public:

   Edisp3(const irfUtil::IrfHdus & edisp_hdus, size_t epoch=0,
          size_t nrow=0);

   Edisp3(const std::string & fitsfile,
          const std::string & extname="ENERGY DISPERSION",
          const std::string & scaling_extname="EDISP_SCALING_PARAMS",
          size_t nrow=0);

   virtual ~Edisp3();

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
      return new Edisp3(*this);
   }

   double scaleFactor(double energy, double costheta) const;

   double evaluate(double emeas, double energy,
                   double theta, double phi, double time, 
                   double * pars) const;

   void renormalize(double logE, double costh, double * params) const {
      (void)(logE);
      (void)(costh);
      (void)(params);
      // Do nothing function to support templated interface expected
      // by EdispInterpolator.  Edisp3 supposedly is automatically
      // correctly normalized.
   }

private:

   std::string m_fitsfile;
   std::string m_extname;
   size_t m_nrow;

   ParTables m_parTables;

   mutable double m_loge_last;
   mutable double m_costh_last;

   mutable double m_pars[10];


   mutable EdispInterpolator * m_interpolator;

   double * pars(double energy, double costh) const;

   std::vector<double> m_scalePars;
   double m_p1;
   double m_p2;
   double m_t0;

   double thibaut_function(double xx, double * pars) const;
   double thibaut_base_function(double xx, double sigma, double kk,
                                double bb, double pp) const;

   void readScaling(const std::string & fitsfile,
                    const std::string & extname);

   class EdispIntegrand {
   public:
      EdispIntegrand(double * pars, double etrue, double scaleFactor,
                     const Edisp3 & self) 
         : m_pars(pars), m_etrue(etrue), m_scaleFactor(scaleFactor),
           m_self(self) {}

      double operator()(double emeas) const {
         double xx((emeas - m_etrue)/m_etrue/m_scaleFactor);
         return m_self.thibaut_function(xx, m_pars)/m_etrue;
      }

   private:
      double * m_pars;
      double m_etrue;
      double m_scaleFactor;
      const Edisp3 & m_self;
   };

};

} // namespace latResponse

#endif // latResponse_Edisp3_h
