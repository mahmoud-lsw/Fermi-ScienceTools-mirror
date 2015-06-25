/** 
 * @file Arg.h
 * @brief Declaration of Arg class
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Arg.h,v 1.1.1.1 2003/08/02 22:14:22 jchiang Exp $
 */

#ifndef optimizers_Arg_h
#define optimizers_Arg_h

namespace optimizers {

/** 
 * @class Arg
 *
 * @brief An abstract class that encapsulates argument type
 * information so that Function's value() and Parameter passing
 * methods can be overloaded transparently.
 *
 * @authors J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Arg.h,v 1.1.1.1 2003/08/02 22:14:22 jchiang Exp $
 */

class Arg {
    
public:

   Arg() {}
   
   virtual ~Arg() {}

protected:

//   Arg() {}

};

} // namespace optimizers

#endif // optimizers_Arg_h
