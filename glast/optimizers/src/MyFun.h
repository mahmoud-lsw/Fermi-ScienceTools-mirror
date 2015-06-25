/** 
 * @file MyFun.h
 * @brief Test function declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/MyFun.h,v 1.4 2015/03/03 18:03:49 jchiang Exp $
 */

#include "optimizers/Function.h"

namespace optimizers {

class Arg;

/** 
 * @class MyFun
 *
 * @brief A simple test function that inherits from Function
 *
 */
    
class MyFun : public Function {

public:

   MyFun();

   ~MyFun(){}

   virtual MyFun * clone() const {
      return new MyFun(*this);
   }

protected:

   virtual double value(Arg &) const;

   virtual double derivByParamImp(Arg & x, const std::string & paramName) const;

};

} // namespace optimizers

