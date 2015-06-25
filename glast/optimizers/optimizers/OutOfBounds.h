/**
 * @file OutOfBounds.h
 * @brief Declaration/definition of OutOfBounds exception class
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/OutOfBounds.h,v 1.2 2003/08/19 17:30:27 cohen Exp $
 */

#ifndef optimizers_OutOfBounds_h
#define optimizers_OutOfBounds_h

#include "optimizers/Exception.h"

namespace optimizers {

/**
 * @class OutOfBounds
 *
 * @brief Exception class used by Parameter objects to help ensure
 * set[True]Value and setBounds methods behave consistently with
 * regard to existing values.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/OutOfBounds.h,v 1.2 2003/08/19 17:30:27 cohen Exp $
 */

class OutOfBounds : public Exception {

public:
   OutOfBounds(const std::string &errorString, double value, 
               double minValue, double maxValue, int code) : 
      Exception(errorString, code), m_value(value), 
      m_minValue(minValue), m_maxValue(maxValue) {}

   ~OutOfBounds() throw() {}
   
   double value() {return m_value;}
   double minValue() {return m_minValue;}
   double maxValue() {return m_maxValue;}
   
   enum ERROR_CODES {VALUE_ERROR, BOUNDS_ERROR};

private:

   double m_value;
   double m_minValue;
   double m_maxValue;

};

} // namespace optimizers

#endif // optimizers_OutOfBounds_h
