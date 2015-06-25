/**
 * @file TRandom4.cxx
 * @brief Use CLHEP random engine instead of ROOT versions.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/src/TRandom4.cxx,v 1.3 2010/06/03 04:13:53 jchiang Exp $
 */

#include "CLHEP/Random/RandFlat.h"

#include "celestialSources/TRandom4.h"

TRandom4::TRandom4(UInt_t seed) {
   SetSeed(seed);
}

TRandom4::~TRandom4() {
}

Double_t TRandom4::Rndm(Int_t i) {
   (void)(i); 
   return CLHEP::RandFlat::shoot();
}

void TRandom4::RndmArray(Int_t n, Float_t *array) {
   for (Int_t i = 0; i < n; i++) {
      array[i] = Rndm(0);
   }
}

void TRandom4::RndmArray(Int_t n, Double_t *array) {
   for (Int_t i = 0; i < n; i++) {
      array[i] = Rndm(0);
   }
}

void TRandom4::SetSeed(UInt_t seed) {
   if (seed != 0) {
      CLHEP::HepRandom hepRandom(seed);
   }
}
