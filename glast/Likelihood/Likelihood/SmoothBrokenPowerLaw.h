/*
 * @file SmoothBrokenPowerLaw.h
 * @brief Declaration for the SmoothBrokenPowerLaw Function class
 * @author Benoit Lott
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/SmoothBrokenPowerLaw.h,v 1.1.1.2.6.3 2015/04/26 06:53:44 jasercio Exp $
 */

#ifndef Likelihood_SmoothBrokenPowerLaw_h
#define Likelihood_SmoothBrokenPowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class SmoothBrokenPowerLaw
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Benoit Lott
   *
   */
  
  class SmoothBrokenPowerLaw : public optimizers::Function {
    
  public:
    
    SmoothBrokenPowerLaw(double Prefactor=10.,
                         double Index1=-2.1, 
                         double Scale=100.,
                         double Index2=-2.1,
                         double BreakValue=1000.,
                         double Beta=0.2);
    
    virtual optimizers::Function * clone() const {
      return new SmoothBrokenPowerLaw(*this);
    }
    
  protected:

     virtual double value(optimizers::Arg&) const;
    
     virtual double derivByParamImp(optimizers::Arg & x,
                                    const std::string & paramName) const;

  };
  
} // namespace Likelihood

#endif // Likelihood_SmoothBrokenPowerLaw_h
