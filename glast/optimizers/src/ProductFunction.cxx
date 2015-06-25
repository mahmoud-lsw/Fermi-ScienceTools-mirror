/** 
 * @file ProductFunction.cxx
 * @brief ProductFunction class implementation
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/ProductFunction.cxx,v 1.2 2015/03/03 18:03:49 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "optimizers/ProductFunction.h"

namespace optimizers {

ProductFunction::ProductFunction(Function & a, Function & b)
   : CompositeFunction(a, b) {
   assert( (a.funcType() == Addend && b.funcType() == Factor) || 
           (a.funcType() == Factor && b.funcType() == Addend) || 
           (a.funcType() == Factor && b.funcType() == Factor) );
   syncParams();
}

void ProductFunction::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                                  bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   std::vector<double> my_derivs;
   if (getFree) {
      m_a->getFreeDerivs(x, my_derivs);
   } else {
      m_a->getDerivs(x, my_derivs);
   }
   for (unsigned int i = 0; i < my_derivs.size(); i++) {
      derivs.push_back(my_derivs[i]*m_b->operator()(x));
   }

   if (getFree) {
      m_b->getFreeDerivs(x, my_derivs);
   } else {
      m_b->getDerivs(x, my_derivs);
   }
   for (unsigned int i = 0; i < my_derivs.size(); i++) {
      derivs.push_back(my_derivs[i]*m_a->operator()(x));
   }
}

} // namespace optimizers
