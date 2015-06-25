/** 
 * @file Psf.cxx
 * @brief Implementation for a test point-spread function class.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/testResponse/src/Psf.cxx,v 1.9 2006/08/14 23:41:22 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Geometry/Vector3D.h"

using CLHEP::RandFlat;
using CLHEP::RandGauss;
using CLHEP::Hep3Vector;

#include "st_facilities/dgaus8.h"
#include "st_facilities/Util.h"

#include "irfInterface/AcceptanceCone.h"

#include "Psf.h"

namespace {
   double erfcc(double x) {
/* (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y.. */
      double t, z, ans;
      
      z=fabs(x);
      t=1.0/(1.0+0.5*z);
      ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
            t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
            t*(-0.82215223+t*0.17087277)))))))));
      return x >= 0.0 ? ans : 2.0-ans;
   }
}

namespace testResponse {

Psf::Psf(const std::vector<double> & psfParams) 
   : m_psfParams(psfParams), m_acceptanceCone(0), m_psfNorm(1.),
     m_haveAngularIntegrals(false) {
   computeCumulativeDist();
}

Psf::~Psf() {
   delete m_acceptanceCone;
}

Psf::Gint Psf::s_gfunc;

Psf::Psf(const Psf &rhs) : IPsf(rhs) {
   m_psfParams = rhs.m_psfParams;
   if (!rhs.m_acceptanceCone) {
      m_acceptanceCone = 0;
   } else {
      m_acceptanceCone 
         = new irfInterface::AcceptanceCone(*rhs.m_acceptanceCone);
   }
   m_psfNorm = rhs.m_psfNorm;
   m_scaledDevs = rhs.m_scaledDevs;
   m_cumDist = rhs.m_cumDist;
   m_psi = rhs.m_psi;
   m_sepMean = rhs.m_sepMean;
   m_angularIntegral = rhs.m_angularIntegral;
   m_needIntegral = rhs.m_needIntegral;
   m_haveAngularIntegrals = rhs.m_haveAngularIntegrals;
}

double Psf::value(const astro::SkyDir &appDir, 
                  double energy, 
                  const astro::SkyDir &srcDir, 
                  const astro::SkyDir &scZAxis,
                  const astro::SkyDir &,
                  double time) const {
   (void)(time);
// Angle between photon and source directions in radians.
   double separation = appDir.difference(srcDir);

// Inclination wrt spacecraft z-axis in radians
   double inc = srcDir.difference(scZAxis);

   return value(separation, energy, inc);
}

astro::SkyDir Psf::appDir(double energy,
                          const astro::SkyDir &srcDir,
                          const astro::SkyDir &scZAxis,
                          const astro::SkyDir &,
                          double time) const {
   (void)(time);
   
   double inclination = srcDir.difference(scZAxis);

   double theta = sepMean(energy, inclination)*drawScaledDev();

   double xi = RandFlat::shoot();
   double phi = 2.*M_PI*xi;

// These might as well have been passed by value...
   astro::SkyDir appDir = srcDir;
   astro::SkyDir zAxis = scZAxis;

// Create an arbitrary x-direction about which to perform the theta
// rotation.
   Hep3Vector srcVec(appDir());
   Hep3Vector arbitraryVec(srcVec.x() + 1., srcVec.y() + 1., srcVec.z() + 1.);
   Hep3Vector xVec = srcVec.cross(arbitraryVec.unit());

   appDir().rotate(theta, xVec).rotate(phi, srcVec);

   return appDir;
}

double Psf::value(double separation, double energy, double theta,
                  double phi, double time) const {
   if (theta < 0) {
      std::ostringstream message;
      message << "testResponse::Psf::value(...):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }
   (void)(phi);
   (void)(time);
   return value(separation*M_PI/180., energy, theta*M_PI/180.);
}

double Psf::value(double separation, double energy, double inc) const {
   return value(separation, sepMean(energy, inc));
}

double Psf::value(double separation, double sepMean) const {
   double gamma = m_psfParams[4];

   double my_value;
   double x = separation/sepMean;
   if (separation != 0) {
      my_value = scaledDist(x)/sepMean;
// Apply normalization and convert to per steradians.
      return my_value/2./M_PI/sin(separation)/m_psfNorm;
   } else {
      my_value = pow((1. + 2.*x/(gamma-1.)), -gamma)/sepMean/sepMean;
      return my_value/2./M_PI/m_psfNorm;
   }
   return 0;
}

double Psf::angularIntegral(double energy,
                            const astro::SkyDir &srcDir,
                            const astro::SkyDir &scZAxis,
                            const astro::SkyDir &,
                            const std::vector<irfInterface::AcceptanceCone *> 
                            &acceptanceCones,
                            double time) {
   double theta = srcDir.difference(scZAxis)*180./M_PI;
   static double phi;
   return angularIntegral(energy, srcDir, theta, phi, acceptanceCones, time);
}

double Psf::angularIntegral(double energy, const astro::SkyDir & srcDir,
                            double theta, double phi, 
                            const std::vector<irfInterface::AcceptanceCone *> 
                            &acceptanceCones, double time) {
   (void)(phi);
   (void)(time);
   if (!m_acceptanceCone || *m_acceptanceCone != *(acceptanceCones[0])) {
      computeAngularIntegrals(acceptanceCones);
      m_haveAngularIntegrals = true;
   }
   double psi = srcDir.difference(m_acceptanceCone->center());
   if (psi >= m_psi.back()) {
      return 0;
   }
   unsigned int ipsi = std::upper_bound(m_psi.begin(), m_psi.end(), psi) 
      - m_psi.begin();
      
   double sep_mean = sepMean(energy, theta*M_PI/180.);
   if (sep_mean < m_sepMean.front() || sep_mean >= m_sepMean.back()) {
      return 0;
   }

   unsigned int isepMean = 
      std::upper_bound(m_sepMean.begin(), m_sepMean.end(), sep_mean) 
      - m_sepMean.begin() - 1;

   unsigned int indx = ipsi*m_sepMean.size() + isepMean;

   if (isepMean == m_sepMean.size() - 1) {
      return psfIntegral(psi, sep_mean);
   }

// Interpolate in sepMean dimension.
   if (m_needIntegral[indx]) {
      m_angularIntegral[indx] = psfIntegral(psi, m_sepMean[isepMean]);
      m_needIntegral[indx] = false;
   }
   if (m_needIntegral[indx+1]) {
      m_angularIntegral[indx+1] = psfIntegral(psi, m_sepMean[isepMean+1]);
      m_needIntegral[indx+1] = false;
   }
   double my_value = (sep_mean - m_sepMean[isepMean])
      /(m_sepMean[isepMean+1] - m_sepMean[isepMean])
      *(m_angularIntegral[indx+1] - m_angularIntegral[indx])
      + m_angularIntegral[indx];
   return my_value;
}

double Psf::angularIntegral(double energy, double theta, 
                            double phi, double radius, double time) const {
   (void)(phi);
   (void)(time);
   double scaledDev = radius*M_PI/180./sepMean(energy, theta*M_PI/180.);
   if (scaledDev >= m_scaledDevs.back()) {
      return 1.;
   } else if (scaledDev <= m_scaledDevs.front()) {
      return 0;
   }
   return st_facilities::Util::interpolate(m_scaledDevs, m_cumDist, scaledDev);
}

double Psf::sepMean(double energy, double inclination) const {
   double mu = cos(inclination);
   double my_value = (m_psfParams[0]*(1. - mu)*(1. - mu) + m_psfParams[1])
      *pow((1./energy + 1./m_psfParams[3]), m_psfParams[2]);
   return my_value;
}

void Psf::computeAngularIntegrals
(const std::vector<irfInterface::AcceptanceCone *> &cones) {
   
// Assume the first acceptance cone is the ROI; ignore the rest.
   if (!m_acceptanceCone) {
      m_acceptanceCone = new irfInterface::AcceptanceCone();
   }
   *m_acceptanceCone = *cones[0];

   unsigned int npsi = 1000;
   unsigned int nsepMean = 200;

   if (!m_haveAngularIntegrals) {
// Set up the array describing the separation between the center of
// the ROI and the source, m_psi
      m_psi.clear();
      double psiMin = 0.;
      double psiMax = 60.*M_PI/180.;    // Is this a large enough angle?
      double psiStep = (psiMax - psiMin)/(npsi - 1.);
      for (unsigned int i = 0; i < npsi; i++) {
         m_psi.push_back(psiStep*i + psiMin);
      }

// The mean separation should be logrithmically spaced because of the
// strong energy dependence of the psf and will have a range given by
// the extremal values of sepMean(energy, inc).
// For now, hard-wire the minimum and maximum energies.
      double emin(20.);
      double emax(2e5);
      m_sepMean.clear();
      double sepMeanMin = sepMean(emax, 0.);
      double sepMeanMax = sepMean(emin, M_PI/2.);
      double sepMeanStep = log(sepMeanMax/sepMeanMin)/(nsepMean - 1.);
      for (unsigned int j = 0; j < nsepMean; j++) {
         m_sepMean.push_back(sepMeanMin*exp(sepMeanStep*j));
      }
   }

// Fill m_angularIntegrals.
   unsigned int npts = m_psi.size()*m_sepMean.size();
   m_angularIntegral.resize(npts);
   m_needIntegral.resize(npts);
   for (unsigned int i = 0; i < npts; i++) {
      m_needIntegral[i] = true;
   }
}

double Psf::psfIntegral(double psi, double sepMean, double roi_radius) {
   if (roi_radius == 0) {
      roi_radius = m_acceptanceCone->radius()*M_PI/180.;
   }
   double one = 1.;
   double mup = cos(roi_radius + psi);
   double mum = cos(roi_radius - psi);
   
   double cp = cos(psi);
   double sp = sin(psi);
   double cr = cos(roi_radius);

   double err = 1e-5;
   long ierr;
   double firstIntegral = 0;

// Check if point is inside or outside the cone            
   if (psi < roi_radius) {
// Integrate the first term...
      s_gfunc = Gint(this, sepMean);
      dgaus8_(&Psf::gfuncIntegrand, &mum, &one, 
              &err, &firstIntegral, &ierr);
   }

// and the second.
   s_gfunc = Gint(this, sepMean, cp, sp, cr);
   double secondIntegral;
   dgaus8_(&Psf::gfuncIntegrand, &mup, &mum, 
           &err, &secondIntegral, &ierr);
   
   double my_integral;
   if (psi < roi_radius) {
      my_integral = firstIntegral + secondIntegral;
   } else {
      my_integral = secondIntegral;
   }
   return my_integral;
}

double Psf::Gint::value(double mu) const {
   double theta = acos(mu);
   double my_value = m_psfObj->value(theta, m_sepMean);

   if (!m_doFirstTerm) {
      double phimin;
      double arg = (m_cr - mu*m_cp)/sqrt(1. - mu*mu)/m_sp;
      if (arg >= 1.) {
         phimin = 0;
      } else if (arg <= -1.) {
         phimin = M_PI;
      } else {
         phimin = acos(arg);
      }
      return 2.*phimin*my_value;
   }
   return 2.*M_PI*my_value;
}

double Psf::scaledDist(double x) const {
   double gamma(m_psfParams[4]);
   return x*pow(1. + 2.*x/(gamma-1.), -gamma);
}

void Psf::computeCumulativeDist() {
   int npts(200);
   double xstep = 6./(npts-1.);
   double xmin(1e-3);
   m_scaledDevs.clear();
   m_scaledDevs.reserve(npts);
   m_scaledDevs.push_back(xmin);
   m_cumDist.clear();
   m_cumDist.reserve(npts);
   m_cumDist.push_back(0);
   for (int i = 1; i < npts; i++) {
      m_scaledDevs.push_back(xmin*pow(10., xstep*i));
      m_cumDist.push_back(m_cumDist[i-1] + (m_scaledDevs[i]-m_scaledDevs[i-1])
                          *(scaledDist(m_scaledDevs[i]) +
                            scaledDist(m_scaledDevs[i-1]))/2.);
   }
   m_psfNorm = m_cumDist.back();
//   std::cerr << "m_psfNorm: " << m_psfNorm << std::endl;
   for (int i = 0; i < npts; i++) {
      m_cumDist[i] /= m_psfNorm;
   }
}

double Psf::drawScaledDev() const {
   double xi = RandFlat::shoot();
   double scaledDev(0);
   try {
      scaledDev 
         = st_facilities::Util::interpolate(m_cumDist, m_scaledDevs, xi);
   } catch (std::exception & eObj) {
      if (st_facilities::Util::
          expectedException(eObj, "abscissa value out-of-range")) {
         scaledDev = 0;
      } else {
         throw;
      }
   }
   return scaledDev;
}

} // namespace testResponse
