/**
 * @file test.cxx
 * @brief Test program for testResponse irfs.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/testResponse/src/test/test.cxx,v 1.6 2007/01/22 00:08:52 jchiang Exp $
 */

#include <cstdlib>

#include <fstream>

#include "astro/SkyDir.h"

#include "irfInterface/IrfsFactory.h"
#include "irfInterface/AcceptanceCone.h"
#include "irfInterface/IrfRegistry.h"

#include "testResponse/IrfLoader.h"

void test_psf(int iargc, char * argv[]);
void test_aeff();
void test_edisp();

int main(int iargc, char * argv[]) {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());
   registry.registerLoader(new testResponse::IrfLoader());
   registry.loadIrfs("TEST");

   test_psf(iargc, argv);
   test_aeff();
   test_edisp();
}

irfInterface::IrfsFactory * irfFactory() {
   return irfInterface::IrfsFactory::instance();
}

void test_edisp() {
   astro::SkyDir scZAxis(180, 0.);
   astro::SkyDir scXAxis(180, 90.);
   astro::SkyDir srcDir(180, 0);
   double energy(100.);

   irfInterface::Irfs * irfs = irfFactory()->create("testIrfs::Front");

   irfInterface::IEdisp * edisp = irfs->edisp();

   int npts(1000);

   std::ofstream outfile("edisp.dat");
   for (int i = 0; i < npts; i++) {
      outfile << edisp->appEnergy(energy, srcDir, scZAxis, scXAxis) 
              << std::endl;
   }
   outfile << edisp->integral(energy*0.8, energy*1.2, energy, 
                              srcDir, scZAxis, scXAxis) << std::endl;
   outfile.close();
}

void test_aeff() {
   irfInterface::Irfs * irfs = irfFactory()->create("testIrfs::Front");
   irfInterface::IAeff * aeff = irfs->aeff();

   double phi(0);
   double theta(0);
   double energy;
   double emin(20);
   double emax(2e5);
   int npts(100);
   double estep = log(emax/emin)/(npts - 1.);
   
   std::ofstream outfile("aeff.dat");
   for (int i = 0; i < npts; i++) {
      energy = emin*exp(estep*i);
      outfile << energy << "  "
              << aeff->value(energy, theta, phi) << std::endl;
   }
   outfile.close();

   std::cout << "Max. effective area: " << aeff->upperLimit() << std::endl;
}

void test_psf(int iargc, char * argv[]) {
   double radius(10.);
   double inclination(0);
   if (iargc >= 2) {
      radius = std::atof(argv[1]);
   }
   if (iargc == 3) {
      inclination = std::atof(argv[2]);
   }

   astro::SkyDir scZAxis(180, 0.);
   astro::SkyDir scXAxis(180, 90.);
   astro::SkyDir srcDir(180+inclination, 0);
   double energy(1000.);

   irfInterface::AcceptanceCone * roiCone = 
      new irfInterface::AcceptanceCone(scZAxis, radius);
   std::vector<irfInterface::AcceptanceCone *> cones;
   cones.push_back(roiCone);

   irfInterface::Irfs * irfs = irfFactory()->create("testIrfs::Front");
   irfInterface::IPsf * psf = irfs->psf();

   int npts(10000);

   int nn(0);
   std::ofstream psf_file("psf.dat");
   for (int i = 0; i < npts; i++) {
      astro::SkyDir appDir = psf->appDir(energy, srcDir, scZAxis, scXAxis);
      double dist = appDir.difference(scZAxis);
      psf_file << dist << "  "
               << appDir.ra() << "  "
               << appDir.dec() << std::endl;
      if (dist <= radius*M_PI/180.) {
         nn++;
      }
   }
   psf_file.close();
   std::cerr << "total within " << radius << " degrees: " << nn << "\n"
             << "integral fraction: " << "\n"
             << psf->angularIntegral(energy, inclination, 0, radius) << "\n"
             << psf->angularIntegral(energy, srcDir, scZAxis, scXAxis,
                                     cones) << "\n"
             << std::endl;
}
