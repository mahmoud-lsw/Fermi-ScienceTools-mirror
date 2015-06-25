/**
 * @file Functor.h
 * @brief Function object interface definition to be used by Amoeba.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Functor.h,v 1.1 2006/07/07 19:41:40 jchiang Exp $
 */

#ifndef optimizers_Functor_h
#define optimizers_Functor_h

#include <vector>

namespace optimizers {

/**
 * @class Functor
 * @brief All function objects to be minimized by Amoeba must derive 
 * from this class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Functor.h,v 1.1 2006/07/07 19:41:40 jchiang Exp $
 */

class Functor {

public:

   Functor() {}

   virtual ~Functor() {}

   virtual double operator()(std::vector<double> & x) = 0;

};

} // namespace optimizers

#endif // optimizers_Functor_h
