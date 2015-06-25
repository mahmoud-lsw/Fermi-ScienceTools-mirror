/**
 * @file Drmnfb.h
 * @brief Declaration for the Drmnfb Optimizer subclass.
 * @author P. Nolan, J. Chiang
 *
 */

#ifndef optimizers_drmnfb_h
#define optimizers_drmnfb_h

#include "optimizers/Optimizer.h"
#include "optimizers/Statistic.h"
#include "optimizers/f2c_types.h"
#include <string>

namespace optimizers {
  
  /** 
   * @class Drmnfb
   *
   * @brief Wrapper class for the MNGB optimization method from 
   * Netlib's PORT library.
   *
   * Finds a local minimum of a continuously differentiable function subject 
   * to simple upper and lower bound constraints. User supplies objective 
   * function. Gradient is calculated by finite differences. Uses a 
   * variant of Newton's method with a quasi-Newton (BFGS) Hessian updating 
   * method, and a model/trust-region technique to aid convergence from 
   * poor starting values.
   *
   * @author P. Nolan, J. Chiang
   * 
   */
  
  /**
   * @file "drmnfb_routines.c"
   *
   * @brief Fortran code for DRMNFB translated by f2c
   *
   The original code for DRMNFB, obtained from Netlib, is in Fortran.
   It was translated using f2c -C++.  The only hand modification required
   was to change \#include "f2c.h" to \#include "f2c/f2c.h" in conformity
   with the way CMT wants files to be laid out.  Most of the code is
   shared with DRMNGB.  Only the unique DRMNFB code is here.
   
  */
  
  class Drmnfb : public Optimizer {
    
  public:
    
    Drmnfb(Statistic &stat) : Optimizer(stat), m_NeedCovariance(true) {}
    
    virtual ~Drmnfb() {}
    
    int find_min(int verbose = 0, double tol = 1e-8, 
		  int tolType = ABSOLUTE);
    int find_min_only(int verbose = 0, double tol = 1e-8, 
		  int tolType = ABSOLUTE);

    //! One-sigma confidence regions based on Hessian, assuming 
    // that this function is a likelihood
    virtual std::vector<double> & getUncertainty(bool useBase=true);

    //! Switch on the covariance matrix calculation. 
    // Save time by not doing this, if covariance is not needed.
    void setNeedCovariance(bool);

    virtual std::ostream& put (std::ostream& s) const;

    enum DrmnfbReturnCodes {
      XCONV = 3, RELCONV, BOTHCONV, ABSCONV, SINGCONV, FALSECONV,
      EVALLIM, ITLIM, STOPX, ALLOCATED = 14, LIVSMALL, LVSMALL,
      RESTARTATT, DNEG, VORANGE, CANTCOMPUTEF = 63, BADPAR, 
      CANTCOMPUTEG, BADARR, BADPAR1, BUGS
    };

  private:
    
    int m_retCode;
    bool m_NeedCovariance;
    int m_evals, m_grads;
    double m_val;

  };

} // namespace optimizers

#endif 


