/**
 * @file IrfLoader.h
 * @brief Concrete derived class of irfInterface::IrfLoader
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/dc1aResponse/dc1aResponse/IrfLoader.h,v 1.2 2007/01/22 00:04:15 jchiang Exp $
 */

#ifndef dc1aResponse_IrfLoader_h
#define dc1aResponse_IrfLoader_h

#include "irfInterface/IrfLoader.h"

namespace dc1aResponse {

/**
 * @class IrfLoader
 * @brief Concrete derived class of irfInterface::IrfLoader
 */

class IrfLoader : public irfInterface::IrfLoader {

public:

   virtual void registerEventClasses() const;

   virtual void loadIrfs() const;

   virtual std::string name() const {
      return "DC1A";
   }

};

} // namespace dc1aResponse

#endif // dc1aResponse_IrfLoader_h
