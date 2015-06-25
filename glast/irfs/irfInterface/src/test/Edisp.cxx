/**
 * @file Edisp.cxx
 * @brief Gaussian energy dispersion function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/src/test/Edisp.cxx,v 1.1 2006/10/23 21:01:06 jchiang Exp $
 */

#include <cmath>

#include "astro/SkyDir.h"

#include "Edisp.h"

Edisp::Edisp(double resolution) 
   : irfInterface::IEdisp(), m_resolution(resolution) {}

double Edisp::value(double appEnergy,
                    double energy, 
                    const astro::SkyDir & srcDir, 
                    const astro::SkyDir & scZAxis,
                    const astro::SkyDir & scXAxis, 
                    double time) const {
   (void)(scXAxis);
   static double phi(0);

   double theta(srcDir.difference(scZAxis)*180./M_PI);

   return value(appEnergy, energy, theta, phi, time);
}

double Edisp::value(double appEnergy, double energy, double theta, 
                    double phi, double time) const {
   (void)(theta);
   (void)(phi);
   (void)(time);

   double sigma(m_resolution*energy);
   double de(appEnergy - energy);

   return std::exp(-de*de/2./sigma/sigma)/std::sqrt(2.*M_PI)/sigma;
}
