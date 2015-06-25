/**
 * @file GaussianSpectrum.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/celestialSources/genericSources/src/GaussianSpectrum.cxx,v 1.1.1.2 2012/09/10 18:39:09 areustle Exp $
 */

#include <cmath>
#include <cstdlib>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/GaussianSpectrum.h"

#include "celestialSources/ConstParMap.h"

ISpectrumFactory & GaussianSpectrumFactory() {
   static SpectrumFactory<GaussianSpectrum> myFactory;
   return myFactory;
}

GaussianSpectrum::GaussianSpectrum(const std::string & paramString) 
   : m_flux(1.), m_mean(1e4), m_sigma(1e3), m_l(0), m_b(0) {

   celestialSources::ConstParMap pars(paramString);
   
   m_flux = pars.value("flux");
   m_mean = pars.value("mean");
   m_sigma = pars.value("sigma");
   m_l = pars.value("glon");
   m_b = pars.value("glat");
}

float GaussianSpectrum::operator()(float xi) const {
   (void)(xi);
   float energy(CLHEP::RandGauss::shoot(m_mean, m_sigma));
   return energy;
}

double GaussianSpectrum::flux(double time) const {
   (void)(time);
   return m_flux;
}

double GaussianSpectrum::solidAngle() const {
   return 1;
}

double GaussianSpectrum::interval(double time) {
   double rate = flux(time)*EventSource::totalArea();
   double xi = CLHEP::RandFlat::shoot();
   return -std::log(1. - xi)/rate;
}

double GaussianSpectrum::energy(double time) {
   (void)(time);
   double xi(0);
   return (*this)(xi);
}

std::pair<double, double> GaussianSpectrum::dir(double energy) {
   (void)(energy);
   return std::make_pair(m_l, m_b);
}
