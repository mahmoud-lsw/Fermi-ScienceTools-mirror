/**
 * @file Registrar.cxx
 * @brief Manage irfInterface::IrfLoader instances for each set of IRFs
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfLoader/src/Registrar.cxx,v 1.14 2015/01/08 23:35:17 jchiang Exp $
 */

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "facilities/commonUtilities.h"

#include "st_facilities/Environment.h"

#include "irfInterface/IrfRegistry.h"

#include "dc1aResponse/IrfLoader.h"
#include "latResponse/IrfLoader.h"
#include "testResponse/IrfLoader.h"

#include "Registrar.h"

namespace irfLoader {

Registrar::Registrar() {
   char * caldb(::getenv("CALDB"));
   /// This calls facilities::commonUtilities::setupEnvironment().
   /// Using this Singleton instance instead of setupEnvironment()
   /// ensures that it will only be called once so that any overridden
   /// values will be preserved.
   st_facilities::Environment::instance();
   if (caldb) {
      /// Use the user-specified environment variable instead of the
      /// one that would otherwise be opaquely forced upon the user by
      /// facilities::commonUtilities::setupEnvironment().
      facilities::commonUtilities::setEnvironment("CALDB", caldb, true);
   }
   
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());

   registry.registerLoader(new dc1aResponse::IrfLoader());
   registry.registerLoader(new latResponse::IrfLoader());
   registry.registerLoader(new testResponse::IrfLoader());

   std::vector<std::string> names(registry.irfNames());
}

} // namespace irfLoader
