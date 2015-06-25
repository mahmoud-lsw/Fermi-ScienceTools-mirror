/** 
 * @file LogNormal.h
 * @brief LogNormal class for spectral modeling
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/LogNormal.h,v 1.1.1.2.6.3 2015/04/26 06:53:44 jasercio Exp $
 */

#ifndef Likelihood_LogNormal_h
#define Likelihood_LogNormal_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/** 
 * @class LogNormal
 *
 */
    
class LogNormal : public optimizers::Function {

public:

   LogNormal(double prefactor=1, double log10_mean=3,
             double log10_sigma=2);

   virtual Function * clone() const {
      return new LogNormal(*this);
   }

protected:

   double value(optimizers::Arg &) const;

   double derivByParamImp(optimizers::Arg & x, 
                          const std::string & paramName) const;

   double integral(optimizers::Arg &, optimizers::Arg &) const {
      return 0;
   }
};

} // namespace Likelihood

#endif // Likelihood_LogNormal_h
