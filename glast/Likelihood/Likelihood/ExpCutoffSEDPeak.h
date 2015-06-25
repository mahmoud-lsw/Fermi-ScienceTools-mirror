/**
 * @file ExpCutoffSEDPeak.cxx
 * @brief Implementation for the ExpCutoff with SED peak energy and flux as variables
 * @author Rolf Buehler
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/ExpCutoffSEDPeak.h,v 1.1.1.1.6.3 2015/04/26 06:53:44 jasercio Exp $
 */

#ifndef Likelihood_ExpCutoffSEDPeak_h
#define Likelihood_ExpCutoffSEDPeak_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class ExpCutoffSEDPeak
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   */
  
  class ExpCutoffSEDPeak : public optimizers::Function {
    
  public:
    
    ExpCutoffSEDPeak(double Fpeak=10.,
                     double Index=-2.1, 
                     double Epeak=1000.);
    
    virtual optimizers::Function * clone() const {
      return new ExpCutoffSEDPeak(*this);
    }
    
  protected:

    double value(optimizers::Arg &) const;
    
    double derivByParamImp(optimizers::Arg & x,
                           const std::string & paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_ExpCutoffSEDPeak_h
