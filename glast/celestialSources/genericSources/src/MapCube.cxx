/**
 * @file MapCube.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/genericSources/src/MapCube.cxx,v 1.15 2012/07/03 20:50:15 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "FitsImage.h"
#include "Util.h"

#include "celestialSources/ConstParMap.h"

#include "genericSources/MapCube.h"

ISpectrumFactory & MapCubeFactory() {
   static SpectrumFactory<MapCube> myFactory;
   return myFactory;
}

MapCube::MapCube(const std::string & paramString) : MapSource() {

   std::string fitsFile;
   bool createSubMap(false);

   if (paramString.find("=") == std::string::npos) {
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      fitsFile = params[1];
      if (params.size() > 2) {
         try {
            m_lonMin = std::atof(params.at(2).c_str());
            m_lonMax = std::atof(params.at(3).c_str());
            m_latMin = std::atof(params.at(4).c_str());
            m_latMax = std::atof(params.at(5).c_str());
            createSubMap = true;
         } catch (...) {
            throw std::runtime_error("Error reading sub-map bounds.\n"
                                     "There must be precisely 4 parameters "
                                     "to describe a sub-map.");
         }
      }
   } else {
      celestialSources::ConstParMap parmap(paramString);
      
      m_flux = parmap.value("flux");
      fitsFile = parmap["fitsFile"];
      if (parmap.size() > 2) {
         try {
            m_lonMin = parmap.value("lonMin");
            m_lonMax = parmap.value("lonMax");
            m_latMin = parmap.value("latMin");
            m_latMax = parmap.value("latMax");
            createSubMap = true;
         } catch (...) {
            throw std::runtime_error("Error reading sub-map bounds.\n"
                                     "They must be specified with names "
                                     "lonMin, lonMax, latMin, latMax.");
         }
      }
   }

   facilities::Util::expandEnvVar(&fitsFile);

   readFitsFile(fitsFile, createSubMap);
   checkForNonPositivePixels(fitsFile);
   readEnergyVector(fitsFile);
   makeCumulativeSpectra();
   std::vector<double> totalCounts(m_solidAngles.size());
   for (unsigned int i = 0; i < m_solidAngles.size(); i++) {
      totalCounts[i] = m_spectra[i].back().first;
   }
   makeIntegralDistribution(totalCounts);

//    std::cerr << "Integral over the map: " 
//              << m_mapIntegral << std::endl;
}

void MapCube::checkForNonPositivePixels(const std::string & fitsFile) const {
   std::vector<double>::const_iterator pixel = m_image.begin();
   for ( ; pixel != m_image.end(); ++pixel) {
      if (*pixel <= 0) {
         std::ostringstream message;
         message << "MapCube: There are negative or zero-valued"
                 << " pixels in the FITS image read from" 
                 << fitsFile;
         throw std::runtime_error(message.str());
      }
   }
}

float MapCube::operator()(float xi) const {
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_integralDist.begin(), m_integralDist.end(), xi);
   unsigned int indx = it - m_integralDist.begin();

   double lon, lat;
   samplePixel(indx, lon, lat);

   if (m_axisTypes[0].find_first_of("R") == 0) { // we have Equatorial coords
      astro::SkyDir myDir(lon, lat);
      m_currentDir = std::make_pair(myDir.l(), myDir.b());
   } else { // assume Galactic coordinates
      m_currentDir = std::make_pair(lon, lat);
   }
   float energy = drawEnergy(m_spectra[indx]);
   return energy;
}

double MapCube::energy(double time) {
   (void)(time);
   double xi = CLHEP::RandFlat::shoot();
   return (*this)(xi);
}

double MapCube::mapValue(unsigned int i, unsigned int j, unsigned int k) {
   unsigned int indx = k*m_lon.size()*m_lat.size() + j*m_lon.size() + i;
   return m_image.at(indx);
}

void MapCube::readEnergyVector(const std::string & fitsFile) {

   std::string routineName("MapCube::readEnergyVector");

   int hdu = genericSources::FitsImage::findHdu(fitsFile, "ENERGIES");
   
   int status(0);
   fitsfile * fptr = 0;

   fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   int hdutype(0);
   fits_movabs_hdu(fptr, hdu, &hdutype, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   genericSources::FitsImage::readColumn(fptr, "Energy", m_energies);

   fits_close_file(fptr, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);
}

void MapCube::makeCumulativeSpectra() {
   m_spectra.clear();
   m_spectra.reserve(m_solidAngles.size());

   for (unsigned int j = 0; j < m_lat.size(); j++) {
      for (unsigned int i = 0; i < m_lon.size(); i++) {
         Spectrum_t counts_spectrum;
         counts_spectrum.reserve(m_energies.size());
         counts_spectrum.push_back(std::make_pair(0, -2));
         for (unsigned int k = 1; k < m_energies.size(); k++) {
            double gamma;
            double counts = counts_spectrum[k-1].first + 
               powerLawIntegral(m_energies.at(k-1), m_energies.at(k),
                                mapValue(i, j, k-1), mapValue(i, j, k),
                                gamma);
            counts_spectrum.push_back(std::make_pair(counts, gamma));
         }
         m_spectra.push_back(counts_spectrum);
      }
   }
}

double MapCube::powerLawIntegral(double x1, double x2,
                                 double y1, double y2, double & gamma) const {
   gamma = std::log(y2/y1)/std::log(x2/x1);
   double n0 = y1/std::pow(x1, gamma);
   double integral;
   if (gamma != 1.) {
      double gp1 = gamma + 1.;
      integral = n0/gp1*(std::pow(x2, gp1) - std::pow(x1, gp1));
   } else {
      integral = n0*std::log(x2/x1);
   }
   return integral;
}

double MapCube::
drawEnergy(const Spectrum_t & spectrum) const {
   double xi = CLHEP::RandFlat::shoot()*spectrum.back().first;
   Spectrum_t::const_iterator it 
      = std::upper_bound(spectrum.begin(), spectrum.end(),
                         std::make_pair(xi, 0), cmpPair);
   int indx = it - spectrum.begin() - 1;
   int nmax = spectrum.size() - 2;
   indx = std::min(std::max(0, indx), nmax);
   double value 
      = genericSources::Util::drawFromPowerLaw(m_energies.at(indx),
                                               m_energies.at(indx+1),
                                               -spectrum.at(indx+1).second);
   return value;
}

bool MapCube::cmpPair(const SpectralPoint_t & x, 
                      const SpectralPoint_t & y) {
   return x.first < y.first;
}
