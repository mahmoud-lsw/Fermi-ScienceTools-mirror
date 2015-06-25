/** 
 * @file Statistic.h
 * @brief Declaration of Statistic base class
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Statistic.h,v 1.5 2015/03/03 18:03:48 jchiang Exp $
 */

#ifndef optimizers_Statistic_h
#define optimizers_Statistic_h

#include <vector>

#include "optimizers/Function.h"

namespace optimizers {

/** 
 * @class Statistic
 *
 * @brief A base class for forcing Function objects to have a
 * convenient interface to use as objective functions, which take no
 * arguments.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Statistic.h,v 1.5 2015/03/03 18:03:48 jchiang Exp $
 */

class Statistic : public Function {
    
public:

   Statistic(const std::string & genericName, unsigned int maxNumParams) 
      : Function(genericName, maxNumParams, "", "dArg", None) {
   }

   virtual ~Statistic() {}

   virtual double value() const = 0;

   virtual void getFreeDerivs(std::vector<double> &derivs) const = 0;

protected:

   Statistic() : Function("Statistic", 0, "", "", None) {}

   Statistic(const Statistic &rhs) : Function (rhs) {}

};

} // namespace optimizers

#endif // optimizers_Statistic_h
