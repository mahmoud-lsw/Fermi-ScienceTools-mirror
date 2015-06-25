/**
 * @file IrfLoader.h
 * @brief Abstract base class for defining the interface expected by
 * IrfRegistry for irfInterface concrete implementations.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/irfInterface/IrfLoader.h,v 1.1 2007/01/22 00:02:58 jchiang Exp $
 */

#ifndef irfInterface_IrfLoader_h
#define irfInterface_IrfLoader_h

#include <string>

namespace irfInterface {

/**
 * @class IrfLoader
 *
 * @brief Abstract base class for defining the interface expected by
 * IrfRegistry for irfInterface concrete implementations.
 */

class IrfLoader {

public:

   IrfLoader() {}

   virtual void registerEventClasses() const = 0;

   virtual void loadIrfs() const = 0;

   virtual std::string name() const = 0;

};

} // namespace irfInterface

#endif // irfInterface_IrfLoader_h
