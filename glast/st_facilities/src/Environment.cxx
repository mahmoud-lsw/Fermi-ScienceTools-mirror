/**
 * @file Environment.cxx
 * @brief Yet another environment interface.  This provides access to
 * facilities::commonUtilities functions that rely on environment
 * variables.  By implementing as a Singleton and providing access to
 * the underlying functions only via the Singleton object, this class
 * ensures that the facilities::commonUtilities::setupEnvironment()
 * function is called without having to burden the clients with this
 * task.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/src/Environment.cxx,v 1.3 2012/11/11 03:37:47 jchiang Exp $
 */

#include "facilities/commonUtilities.h"

#include "st_facilities/Environment.h"

namespace st_facilities {

Environment * Environment::s_instance(0);

Environment & Environment::instance() {
   if (s_instance == 0) {
      s_instance = new Environment();
   }
   return *s_instance;
}

Environment::Environment() {
   facilities::commonUtilities::setupEnvironment();
}

std::string Environment::dataPath(const std::string & package) {
   instance();
   return facilities::commonUtilities::getDataPath(package);
}

std::string Environment::getEnv(const std::string & envvar) {
   instance();
   return facilities::commonUtilities::getEnvironment(envvar);
}

std::string Environment::packagePath(const std::string & package) {
   instance();
   return facilities::commonUtilities::getPackagePath(package);
}

std::string Environment::xmlPath(const std::string & package) {
   instance();
   return facilities::commonUtilities::getXmlPath(package);
}

} // namespace st_facilities
