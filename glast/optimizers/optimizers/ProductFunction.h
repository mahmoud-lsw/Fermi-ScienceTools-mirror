/** 
 * @file ProductFunction.h
 * @brief Declaration of ProductFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/ProductFunction.h,v 1.2 2015/03/03 18:03:48 jchiang Exp $
 */

#ifndef optimizers_ProductFunction_h
#define optimizers_ProductFunction_h

#include "optimizers/CompositeFunction.h"

namespace optimizers {
/** 
 * @class ProductFunction
 *
 * @brief A Function that returns the product of two Functions
 *
 */
    
class ProductFunction : public CompositeFunction {
public:

   ProductFunction(Function & a, Function & b);

   virtual Function * clone() const {
      return new ProductFunction(*this);
   }

protected:

   double value(Arg &x) const {
      return m_a->operator()(x)*m_b->operator()(x);
   }

   void fetchDerivs(Arg &x, std::vector<double> &derivs, bool getFree) const;

};

} // namespace optimizers

#endif // optimizers_ProductFunction_h
