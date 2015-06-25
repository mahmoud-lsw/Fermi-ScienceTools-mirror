/**
 * @file IrfLoader.cxx
 * @brief Concrete implementation of irfInterface/IrfLoader
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/testResponse/src/IrfLoader.cxx,v 1.5 2007/01/30 03:53:26 jchiang Exp $
 */

#include "irfInterface/IrfRegistry.h"

#include "testResponse/IrfLoader.h"

namespace testResponse {

void load_irfs();

void IrfLoader::registerEventClasses() const {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());
   const char * class_names[] = {"testIrfs::Front", "testIrfs::Back"};
   std::vector<std::string> classNames(class_names, class_names + 2);
   registry.registerEventClasses("TEST", classNames);
   registry.registerEventClass("TESTF", "testIrfs::Front");
   registry.registerEventClass("TESTB", "testIrfs::Back");
}

void IrfLoader::loadIrfs() const {
   load_irfs();
}

} // namespace testResponse
