/**
 * @file main.cxx
 * @brief Test program to exercise observationSim interface as a
 * prelude to the O2 tool.
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/src/test/main.cxx,v 1.43 2013/02/11 17:10:16 jchiang Exp $
 */
#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cstdlib>

#include "facilities/commonUtilities.h"

#include "astro/SkyDir.h"

#include "st_facilities/Environment.h"

#include "irfInterface/IrfsFactory.h"
#include "irfLoader/Loader.h"

#include "celestialSources/SpectrumFactoryLoader.h"

#include "dataSubselector/Cuts.h"

#include "observationSim/Simulator.h"
#include "observationSim/EventContainer.h"
#include "observationSim/ScDataContainer.h"
#include "LatSc.h"

void help();

void load_sources();

int main(int iargc, char * argv[]) {
#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

   try {
// Create list of xml input files for source definitions.
   std::vector<std::string> fileList;
   std::string xml_list(facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("observationSim"), "obsSim_source_library.xml"));
   fileList.push_back(xml_list);
   xml_list = facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("observationSim"), "3EG_catalog_20-1e6MeV.xml");
   fileList.push_back(xml_list);
#ifndef BUILD_WITHOUT_ROOT
   xml_list = facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("GRB"), "GRB_user_library.xml");
   fileList.push_back(xml_list);
#endif
   xml_list = facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("observationSim"), "time_source.xml");
   fileList.push_back(xml_list);

   load_sources();

// Parse the command line arguments.
//
// The first argument will be the total number of photons or the total
// simulation time in seconds.
//
   long count;
   if (iargc > 1) {
      count = static_cast<long>(std::atof(argv[1]));
   } else {
      count = 1000;
   }

// All subsequent arguments are either option flags or the names of
// sources.
//
   bool useSimTime(false);
   bool useCombined(true);
   std::vector<std::string> sourceNames;
   if (iargc > 2) {
      for (int i = 2; i < iargc; i++) {
         std::string argString = argv[i];
         if (argString == "-t") {
// Interpret arg[1] as elapsed time in seconds.
            useSimTime = true;
         } else if (argString == "-fb") {
// Use Front/Back responses instead of Combined.
            useCombined = false;
         } else {
// Assume the next argument is a source name or a request for help.
            if (argString == "help") {
               help();
               return 0;
            }
            sourceNames.push_back(argString);
         }
      }
   } else {
      sourceNames.push_back("all_3EG_sources");
//      sourceNames.push_back("One_GRB_each_10Minutes");
//      sourceNames.push_back("galdiffusemap");
   }

// Create the Simulator object
   observationSim::Simulator my_simulator(sourceNames, fileList, 1.21);

// Allow for multiple IRFs.
   irfLoader::Loader::go();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
   std::vector<irfInterface::Irfs *> respPtrs;

   respPtrs.push_back(myFactory->create("DC1A::Front"));
   respPtrs.push_back(myFactory->create("DC1A::Back"));

   dataSubselector::Cuts * cuts(new dataSubselector::Cuts);
   cuts->setIrfs("DC1A");

// Generate the events and spacecraft data.
   observationSim::EventContainer events("test_events", "EVENTS", cuts);
   observationSim::ScDataContainer scData("test_scData", "SC_DATA");

// The spacecraft object.
   observationSim::Spacecraft *spacecraft = new observationSim::LatSc();

// Use simulation time rather than total counts if desired.
   if (useSimTime) {
      std::cout << "Generating events for a simulation time of "
                << count << " seconds....";
      my_simulator.generateEvents(static_cast<double>(count), events, 
                                  scData, respPtrs, spacecraft);
   } else {
      std::cout << "Generating " << count << " events....";
      my_simulator.generateEvents(count, events, scData, respPtrs, spacecraft);
   }
   std::cout << "Done." << std::endl;
   } catch (std::exception & eObj) {
      std::cout << eObj.what() << std::endl;
      return 1;
   }
}

void help() {
   std::cerr << "usage: <program name> <counts> [<options> <sourceNames>]\n"
             << "options: \n"
             << "  -t interpret counts as elapsed time in seconds\n" 
             << "  -fb use Front/Back IRFs\n"
             << std::endl;
}

void load_sources() {
   SpectrumFactoryLoader foo;
}
