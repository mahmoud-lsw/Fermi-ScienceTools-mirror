/**
 * @file test.cxx
 * @brief Code to check SpectrumFactoryLoader class
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/src/test/test.cxx,v 1.1 2005/04/26 23:51:15 jchiang Exp $
 */

#include <iostream>
#include <string>
#include <vector>

#include "celestialSources/SpectrumFactoryLoader.h"

int main() {
   
   SpectrumFactoryLoader loader;

   const std::vector<std::string> & names = loader.names();
   std::cout << "Loaded the following factories: \n";
   for (unsigned int i = 0; i < names.size(); i++) {
      std::cout << names.at(i) << "\n";
   }

}
