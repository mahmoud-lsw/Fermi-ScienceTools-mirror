/** 
 * @file SrcArg.h
 * @brief Declaration of SrcArg class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SrcArg.h,v 1.7 2003/08/06 20:52:03 jchiang Exp $
 */

#ifndef Likelihood_SrcArg_h
#define Likelihood_SrcArg_h

#include "optimizers/Arg.h"
#include "Likelihood/Source.h"

namespace Likelihood {

/** 
 * @class SrcArg
 *
 * @brief Concrete Arg subclass for encapsulating data of type Source.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SrcArg.h,v 1.7 2003/08/06 20:52:03 jchiang Exp $
 */

class SrcArg : public optimizers::Arg {
    
public:
   
   SrcArg(Source *src) : m_val(src) {}
   virtual ~SrcArg() {}

   Source *getValue() const {return m_val;}

private:

   Source *m_val;

};

} // namespace Likelihood

#endif // Likelihood_SrcArg_h
