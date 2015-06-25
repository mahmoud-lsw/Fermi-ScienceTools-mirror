/** 
 * @file BrokenPowerLaw.h
 * @brief Declaration for the BrokenPowerLaw Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/BrokenPowerLaw.h,v 1.4 2015/03/03 18:03:49 jchiang Exp $
 */

#ifndef optimizers_BrokenPowerLaw_h
#define optimizers_BrokenPowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace optimizers {

/** 
 * @class BrokenPowerLaw
 *
 * @brief A broken power-law function
 *
 */
    
class BrokenPowerLaw : public Function {

public:

   BrokenPowerLaw(double Prefactor=1., double Index1=-1.5,
                  double Index2=-2.5, double BreakValue=1000.);

   virtual Function * clone() const {
      return new BrokenPowerLaw(*this);
   }

protected:

   double value(Arg & xarg) const;

   double derivByParamImp(Arg & xarg, const std::string & paramName) const;

};

} // namespace optimizers

#endif // optimizers_BrokenPowerLaw_h
