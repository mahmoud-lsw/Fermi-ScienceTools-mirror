/**
 * @file IrfLoader.h
 * @brief Concrete derived class of irfInterface::IrfLoader
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/testResponse/testResponse/IrfLoader.h,v 1.2 2007/01/22 00:08:52 jchiang Exp $
 */

#ifndef testResponse_IrfLoader_h
#define testResponse_IrfLoader_h

#include "irfInterface/IrfLoader.h"

namespace testResponse {

/**
 * @class IrfLoader
 * @brief Concrete derived class of irfInterface::IrfLoader
 */

class IrfLoader : public irfInterface::IrfLoader {

public:

   virtual void registerEventClasses() const;

   virtual void loadIrfs() const;

   virtual std::string name() const {
      return "TEST";
   }

};

} // namespace testResponse

#endif // testResponse_IrfLoader_h
