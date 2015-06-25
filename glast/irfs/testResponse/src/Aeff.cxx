/** 
 * @file Aeff.cxx
 * @brief Implementation for the effective Area class.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/testResponse/src/Aeff.cxx,v 1.5 2006/08/14 23:41:21 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "Aeff.h"

namespace testResponse {

double Aeff::value(double energy, 
                   const astro::SkyDir &srcDir, 
                   const astro::SkyDir &scZAxis,
                   const astro::SkyDir &scXAxis, 
                   double time) const {
   (void)(scXAxis);
   double theta = srcDir.difference(scZAxis)*180./M_PI;
   return value(energy, theta, 0, time);
}

double Aeff::value(double energy, double theta, double phi,
                   double time) const {
   (void)(phi);
   (void)(time);
   if (theta < 0) {
      std::ostringstream message;
      message << "testResponse::Aeff::value(double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double prefactor = 1.25*(cos(theta*M_PI/180.) - 0.2);
   if (prefactor > 0) {
      double aeff = m_params[0]
         /(1. + exp(-m_params[2]*log(energy/m_params[1])));
      return prefactor*aeff;
   }
   return 0;
}

double Aeff::upperLimit() const {
   return m_params[0];
}

} // namespace testResponse
