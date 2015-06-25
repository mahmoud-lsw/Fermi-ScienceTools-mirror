/**
 * @file Aeff.h
 * @brief Aeff class declaration for post-handoff review IRFs.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/Aeff.h,v 1.1.1.6.2.4 2015/04/26 06:53:49 jasercio Exp $
 */
  
#ifndef latResponse_Aeff_h
#define latResponse_Aeff_h

#include <string>
#include <utility>
#include <vector>

#include "irfInterface/IAeff.h"

#include "latResponse/FitsTable.h"

namespace irfUtil {
   class IrfHdus;
}

namespace latResponse {

   class ParTables;

/**
 * @class Aeff
 * @brief Effective area class declaration.
 */
  
class Aeff : public irfInterface::IAeff {

public:

   Aeff(const irfUtil::IrfHdus & irf_hdus, size_t iepoch=0, 
        size_t nrow=0);

   Aeff(const std::string & fitsfile, 
        const std::string & extname="EFFECTIVE AREA",
        size_t nrow=0);

   Aeff(const Aeff & other);
   
   Aeff & operator=(const Aeff & rhs);

   virtual ~Aeff() {}
   
   virtual double value(double energy, 
                        const astro::SkyDir & srcDir, 
                        const astro::SkyDir & scZAxis,
                        const astro::SkyDir & scXAxis,
                        double time=0) const;
   
   virtual double value(double energy, double theta, double phi,
                        double time=0) const;
   
   virtual irfInterface::IAeff * clone() {
      return new Aeff(*this);
   }

   virtual double upperLimit() const;

   double max_phi_modulation() const;

   std::pair<double, double> pars(double logE, double costh,
                                  bool interpolate=false) const;

   double phi_modulation(double logE, double costheta, double phi,
                         bool interpolate) const;

private:

   FitsTable m_aeffTable;

   ParTables * m_phiDepPars;

   double phi_modulation(double par0, double par1, double phi) const;

};

} // namespace latResponse

#endif // latResponse_Aeff_h
