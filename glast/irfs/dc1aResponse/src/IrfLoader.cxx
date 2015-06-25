/**
 * @file IrfLoader.cxx
 * @brief Concrete implementation of irfInterface/IrfLoader
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/dc1aResponse/src/IrfLoader.cxx,v 1.3 2007/01/23 07:12:16 jchiang Exp $
 */

#include "irfInterface/IrfRegistry.h"

#include "dc1aResponse/IrfLoader.h"

namespace dc1aResponse {

void load_irfs();

void IrfLoader::registerEventClasses() const {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());
   const char * class_names[] = {"DC1A::Front", "DC1A::Back"};
   std::vector<std::string> classNames(class_names, class_names + 2);
   registry.registerEventClasses("DC1A", classNames);
   registry.registerEventClass("DC1AF", "DC1A::Front");
   registry.registerEventClass("DC1AB", "DC1A::Back");
}

void IrfLoader::loadIrfs() const {
   load_irfs();
}

} // namespace dc1aResponse
