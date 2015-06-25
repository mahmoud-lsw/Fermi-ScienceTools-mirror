/**
 * @file Amoeba.h
 * @brief Class to minimize a functor using Nelder-Mead.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Amoeba.h,v 1.6 2009/03/23 23:43:44 jchiang Exp $
 */

#ifndef optimizers_Amoeba_h
#define optimizers_Amoeba_h

#include "optimizers/Functor.h"

namespace optimizers {

/**
 * @class Amoeba
 *
 * @brief This class minimizes a function object that takes
 * std::vector<double> of parameter values using the Nelder-Mead
 * algorithm implemented in NR.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Amoeba.h,v 1.6 2009/03/23 23:43:44 jchiang Exp $
 */

class Amoeba {

public:

   Amoeba(Functor & functor, const std::vector<double> & params,
          double step=0.1, bool addstep=false) 
      : m_functor(functor), m_npars(params.size()) {
      buildSimplex(params, step, addstep);
   }

   ~Amoeba() {}

   /// @return The estimate of the minimum value of the functor.
   /// @param params The parameter values at the minimum estimate.
   /// @param tol The absolute tolerance for convergence.
   double findMin(std::vector<double> & params, double tol=1e-2,
                  bool abstol=true);

private:

   Functor & m_functor;
   size_t m_npars;
   std::vector< std::vector<double> > m_simplex;

   /// @param params An ndim vector of parameter values as a candidate
   /// starting point for the minimization.
   /// @param frac This sets the size of the initial simplex.
   void buildSimplex(const std::vector<double> & params,
                     double frac=0.1, bool addfrac=false);

};

} // namespace optimizers

#endif // optimizers_Amoeba_h
