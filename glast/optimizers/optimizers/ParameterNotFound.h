/**
 * @file ParameterNotFound.h
 * @brief Declaration and definition of ParameterNotFound exception class.
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/ParameterNotFound.h,v 1.1.1.1 2003/08/02 22:14:23 jchiang Exp $
 */

#ifndef optimizers_ParameterNotFound_h
#define optimizers_ParameterNotFound_h

#include <sstream>
#include "optimizers/Exception.h"

namespace optimizers {

/**
 * @class ParameterNotFound
 *
 * @brief A class that returns a standard error message for
 * Parameters looked for but not found in the desired Function.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/ParameterNotFound.h,v 1.1.1.1 2003/08/02 22:14:23 jchiang Exp $
 */

class ParameterNotFound : public Exception {

public:
   ParameterNotFound(const std::string &paramName, 
                     const std::string &funcName,
                     const std::string &routineName) {
      std::ostringstream errorMessage;
      errorMessage << "Function::" << routineName << ": \n"
                   << "A Parameter named " << paramName
                   << " is not a Parameter of Function "
                   << funcName << "\n";
      m_what = errorMessage.str();
      m_code = 0;
   }
};

} // namespace optimizers

#endif // optimizers_ParameterNotFound_h
