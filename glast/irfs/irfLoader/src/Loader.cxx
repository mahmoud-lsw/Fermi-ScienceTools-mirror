/**
 * @file Loader.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfLoader/src/Loader.cxx,v 1.24 2007/06/16 23:11:04 jchiang Exp $
 */

#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "irfInterface/IrfRegistry.h"

#define ST_DLL_EXPORTS
#include "irfLoader/Loader.h"
#undef ST_DLL_EXPORTS

#include "Registrar.h"

namespace irfLoader {

Registrar * Loader::s_registrar(new Registrar());

void Loader::go(const std::string & name) {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());
   std::vector<std::string> irfNames(registry.irfNames());
   if (std::find(irfNames.begin(), irfNames.end(), name) != irfNames.end()) {
      if (name != "EGRET" || ::getenv("CALIB_DIR")) {
         registry.loadIrfs(name);
      }
   } else {
      throw std::invalid_argument("Invalid IRF named " + name + ".");
   }
}

void Loader::go(const std::vector<std::string> & irfsNames) {
   for (std::vector<std::string>::const_iterator name(irfsNames.begin());
        name != irfsNames.end(); ++name) {
      go(*name);
   }
}

void Loader::go() {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());
   std::vector<std::string> names(registry.irfNames());
   for (size_t i(0); i < names.size(); i++) {
      go(names.at(i));
   }
}

const std::map<std::string, std::vector<std::string> > & Loader::respIds() {
   return irfInterface::IrfRegistry::instance().respIds();
}

void Loader_go() {
   Loader::go();
}

const std::map<std::string, std::vector<std::string> > & Loader_respIds() {
   return Loader::respIds();
}

} // namespace irfLoader
