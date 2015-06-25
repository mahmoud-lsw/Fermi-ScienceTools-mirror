// -*- C++ -*-
// $Id: KroneckerDelta.cc,v 1.1.1.1 2012/02/21 18:10:26 elwinter Exp $
#include "CLHEP/GenericFunctions/KroneckerDelta.hh"
#include "CLHEP/GenericFunctions/FixedConstant.hh"
#include <assert.h>

namespace Genfun {
FUNCTION_OBJECT_IMP(KroneckerDelta)

KroneckerDelta::KroneckerDelta()
{}

KroneckerDelta::KroneckerDelta(const KroneckerDelta & right) 
{
}

KroneckerDelta::~KroneckerDelta() {
}

double KroneckerDelta::operator() (double x) const {
  if (x < 0.5 && x > -0.5) return 1;
  return 0;
}


Derivative KroneckerDelta::partial(unsigned int index) const {
  assert(index==0);
  const AbsFunction & fPrime  = FixedConstant(0);;
  return Derivative(&fPrime);
}

} // namespace Genfun
