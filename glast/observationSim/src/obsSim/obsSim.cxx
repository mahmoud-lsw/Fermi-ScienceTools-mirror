/**
 * @file obsSim.cxx
 * @brief Observation simulator using instrument response functions.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/observationSim/src/obsSim/obsSim.cxx,v 1.1.1.22.2.4 2015/04/26 06:54:00 jasercio Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cstdlib>
#include <memory>
#include <stdexcept>

#include "CLHEP/Random/Random.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "facilities/commonUtilities.h"
#include "facilities/Timestamp.h"
#include "facilities/Util.h"

#include "astro/GPS.h"
#include "astro/JulianDate.h"
#include "astro/SkyDir.h"

#include "irfInterface/IrfsFactory.h"
#include "irfUtil/Util.h"
#include "irfLoader/Loader.h"

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "flux/Spectrum.h"

#include "dataSubselector/Cuts.h"

#include "celestialSources/SpectrumFactoryLoader.h"

#include "observationSim/Simulator.h"
#include "observationSim/EventContainer.h"
#include "observationSim/ScDataContainer.h"

#include "LatSc.h"

using st_facilities::Util;

class ObsSim : public st_app::StApp {
public:
   ObsSim() : st_app::StApp(), m_pars(st_app::StApp::getParGroup("gtobssim")),
              m_simulator(0), 
              m_formatter(new st_stream::StreamFormatter("gtobssim", "", 2)) {
      setVersion(s_cvs_id);
   }
   virtual ~ObsSim() throw() {
      try {
         delete m_simulator;
         delete m_formatter;
         for (unsigned int i = 0; i < m_respPtrs.size(); i++) {
            delete m_respPtrs[i];
         }
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) { 
      }
   }
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   double m_count;
   std::vector<std::string> m_xmlSourceFiles;
   std::vector<std::string> m_srcNames;
   std::vector<irfInterface::Irfs *> m_respPtrs;
   observationSim::Simulator * m_simulator;
   st_stream::StreamFormatter * m_formatter;
   double m_tstart;

   void promptForParameters();
   void checkOutputFiles();
   void setRandomSeed();
   void createFactories();
   void setXmlFiles();
   void readSrcNames();
   void createResponseFuncs();
   void createSimulator();
   void generateData();
   void saveEventIds(const observationSim::EventContainer & events) const;
   double maxEffArea() const;
   void get_tstart(std::string scfile, const std::string & sctable);

   static std::string s_cvs_id;
};

st_app::StAppFactory<ObsSim> myAppFactory("gtobssim");

std::string ObsSim::s_cvs_id("$Name: ScienceTools-v10r0p5-fssc-20150518 $");

void ObsSim::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ObsSim::run() {
   promptForParameters();
   checkOutputFiles();
   setRandomSeed();
   createFactories();
   setXmlFiles();
   readSrcNames();
   createResponseFuncs();
   createSimulator();
   generateData();
   m_formatter->info() << "Done." << std::endl;
}

void ObsSim::promptForParameters() {
   m_pars.Prompt("infile");
   m_pars.Prompt("srclist");
   m_pars.Prompt("scfile");
   m_pars.Prompt("evroot");
   m_pars.Prompt("simtime");
   m_pars.Prompt("tstart");
   m_pars.Prompt("use_ac");
   if (m_pars["use_ac"]) {
      m_pars.Prompt("ra");
      m_pars.Prompt("dec");
      m_pars.Prompt("radius");
   }
   m_pars.Prompt("irfs");
   m_pars.Prompt("seed");
   m_pars.Save();
   m_count = m_pars["simtime"];
}

void ObsSim::checkOutputFiles() {
   bool clobber = m_pars["clobber"];
   if (!clobber) {
      std::string prefix = m_pars["evroot"];
      std::string file = prefix + "_events_0000.fits";
      if (st_facilities::Util::fileExists(file)) {
         m_formatter->err() << "Output file " << file  << " already exists,\n"
                            << "and you have set 'clobber' to 'no'.\n"
                            << "Please provide a different output file prefix."
                            << std::endl;
         std::exit(1);
      }
      file = prefix + "_scData_0000.fits";
      if (st_facilities::Util::fileExists(file)) {
         m_formatter->err() << "Output file " << file << " already exists,\n"
                            << "and you have set 'clobber' to 'no'.\n"
                            << "Please provide a different output file prefix."
                            << std::endl;
         std::exit(1);
      }
   }
}

void ObsSim::setRandomSeed() {
// Set the random number seed in the CLHEP random number engine.
// We only do this once per run, so we set it using the constructor.
// See <a href="http://wwwasd.web.cern.ch/wwwasd/lhc++/clhep/doxygen/html/Random_8h-source.html">CLHEP/Random/Random.h</a>.
   CLHEP::HepRandom hepRandom(m_pars["seed"]);
}

void ObsSim::createFactories() {
   SpectrumFactoryLoader foo;
}

void ObsSim::setXmlFiles() {
   m_xmlSourceFiles.clear();
// observationSim::Simulator requires a specific "TimeCandle" source,
// so time_source.xml must always be loaded.
   m_xmlSourceFiles.push_back(facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("observationSim"), "time_source.xml"));

// Fetch any user-specified xml file of flux-style source definitions,
// replacing the default list.
   std::string xmlFiles = m_pars["infile"];
   if (xmlFiles == "none" || xmlFiles == "") { // use the default
      xmlFiles = facilities::commonUtilities::joinPath(st_facilities::Environment::xmlPath("observationSim"), "xmlFiles.dat");
   }
   facilities::Util::expandEnvVar(&xmlFiles);
   if (Util::fileExists(xmlFiles)) {
      if (Util::isXmlFile(xmlFiles)) {
         m_xmlSourceFiles.push_back(xmlFiles);
      } else {
         std::vector<std::string> files;
         Util::readLines(xmlFiles, files, "#", true);
         for (unsigned int i=0; i < files.size(); i++) {
            facilities::Util::expandEnvVar(&files[i]);
            if (Util::fileExists(files[i])) {
               m_xmlSourceFiles.push_back(files[i]);
            } else {
               m_formatter->info() << "File not found: " 
                                   << files[i] << std::endl;
            }
         }
      } 
   } else {
      throw std::invalid_argument("List of XML files not found: " + xmlFiles);
   }
}

void ObsSim::readSrcNames() {
   std::string srcListFile = m_pars["srclist"];
   if (Util::fileExists(srcListFile)) { 
      Util::readLines(srcListFile, m_srcNames, "#", true);
      if (m_srcNames.size() == 0) {
         throw std::invalid_argument("No sources given in " + srcListFile);
      }
   } else {
      throw std::invalid_argument("Source_list file " + srcListFile
                                  + " doesn't exist.");
   }
}   

void ObsSim::createResponseFuncs() {
   irfLoader::Loader::go();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();

   std::string responseFuncs = m_pars["irfs"];
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(responseFuncs, "_", tokens);
   if (tokens.size() > 3) {
      throw std::runtime_error("Invalid IRF designation: " + responseFuncs);
   }
   std::string evtype = m_pars["evtype"];
   if (evtype != "" && evtype !="none") {
      responseFuncs += (" (" + evtype + ")");
   }
   m_formatter->info(3) << "Using irfs: " << responseFuncs << std::endl;

   if (responseFuncs == "none") {
      m_respPtrs.clear();
      return;
   }

   typedef std::map< std::string, std::vector<std::string> > respMap;
   const respMap & responseIds = irfLoader::Loader::respIds();
   respMap::const_iterator it;
   if ( (it = responseIds.find(responseFuncs)) != responseIds.end() ) {
      const std::vector<std::string> & resps = it->second;
      for (size_t i = 0; i < resps.size(); i++) {
         irfInterface::Irfs * irf(myFactory->create(resps[i]));
         m_formatter->info(3) << "Adding IRF: " << resps[i];
         m_formatter->info(4) << ", with event_type bit " << irf->irfID();
         m_formatter->info(3) << std::endl;
         m_respPtrs.push_back(irf);
      }
   } else {
      std::ostringstream message;
      message << "Invalid response function choice: " << responseFuncs << "\n"
              << "Valid choices are \n";
      for (it = responseIds.begin(); it != responseIds.end(); ++it) {
         message << "   " << it->first << "\n";
      }
      throw std::invalid_argument(message.str());
   }
}   

void ObsSim::createSimulator() {
   double totalArea(maxEffArea());
//   std::cout << "total area: " << totalArea << std::endl;
   std::string pointingHistory = m_pars["scfile"];
   std::string sctable = m_pars["sctable"];
   try {
      double value = m_pars["tstart"];
      m_tstart = value;
   } catch (...) { // Assume tstart = INDEF
      if (pointingHistory == "none") {
         m_tstart = 0;
      } else { // Use TSTART from scfile
         get_tstart(pointingHistory, sctable);
      }
   }
   std::string startDate = m_pars["startdate"];
   facilities::Timestamp start(startDate);
   double offset((astro::JulianDate(start.getJulian()) 
                  - astro::JulianDate::missionStart())
                 *astro::JulianDate::secondsPerDay);
   m_tstart += offset;
   Spectrum::setStartTime(offset);
   double maxSimTime = 3.155e8;
   try {
      maxSimTime = m_pars["maxtime"];
   } catch (std::exception &) {
   }
   m_simulator = new observationSim::Simulator(m_srcNames, m_xmlSourceFiles, 
                                               totalArea, m_tstart,
                                               pointingHistory, maxSimTime,
                                               offset);

   int id_offset = m_pars["offset"];
   m_simulator->setIdOffset(id_offset);

   if (pointingHistory == "none" || pointingHistory == "") {
      try {
         double rocking_angle = m_pars["rockangle"];
         m_simulator->setRocking(3, rocking_angle);
      } catch (...) { // rockangle = INDEF
         // do nothing (i.e., leave at default rocking of 35 deg)
      }
   }
}

void ObsSim::generateData() {
   long nMaxRows = m_pars["maxrows"];
   std::string prefix = m_pars["evroot"];
   std::string ev_table = m_pars["evtable"];
   dataSubselector::Cuts * cuts = new dataSubselector::Cuts;
   cuts->addRangeCut("ENERGY", "MeV", m_pars["emin"], m_pars["emax"]);

   // Setting the irfs also sets the cut on CONVERSION_TYPE and the 
   // bit that is set in the EVENT_CLASS variable.
   std::string irfs = m_pars["irfs"];
   cuts->setIrfs(irfs);
   // Remove the VersionCut containing the IRF_VERSION since that should
   // not appear in an FT1 file.
   cuts->removeVersionCut("IRF_VERSION");

   // If evtype is set and an EVENT_TYPE bit-mask cut is not present,
   // add such a cut corresponding to evtype.
   dataSubselector::BitMaskCut * evtype_cut(0);
   try {
      std::string evtype = m_pars["evtype"];
      if (evtype == "none") {
         throw hoops::Hexception(12, "Don't add this cut for Front/Back", 
                                 "", 0);
      }
      if (!(evtype_cut = cuts->bitMaskCut("EVENT_TYPE"))) {
         // Get the inverse mapping
         std::map<std::string, std::pair<unsigned int, std::string> > 
            event_type_mapping;
         std::vector<std::string> partitions;
         std::map<std::string, unsigned int> bitmask_by_partition;
         irfUtil::Util::get_event_type_mapping(irfs,
                                               event_type_mapping,
                                               partitions,
                                               bitmask_by_partition);
         cuts->addBitMaskCut("EVENT_TYPE", bitmask_by_partition[evtype],
                             cuts->pass_ver());
      }
      delete evtype_cut;
   } catch (hoops::Hexception &) {
      // evtype not set so do nothing
   }

   if (m_pars["use_ac"]) {
      cuts->addSkyConeCut(m_pars["ra"], m_pars["dec"], m_pars["radius"]);
   }
   double start_time(m_tstart);
   double stop_time;
   if (m_pars["nevents"]) {
      stop_time = start_time;
   } else {
      double sim_time(m_pars["simtime"]);  // yes, this is BS.
      stop_time = start_time + sim_time;
   }
   bool applyEdisp = m_pars["edisp"];
   observationSim::EventContainer events(prefix + "_events", ev_table,
                                         cuts, nMaxRows,
                                         start_time, stop_time, applyEdisp,
                                         &m_pars);
   events.setAppName("gtobssim");
   events.setVersion(getVersion());
   std::string pointingHistory = m_pars["scfile"];
   facilities::Util::expandEnvVar(&pointingHistory);
   bool writeScData = (pointingHistory == "" || pointingHistory == "none"
                       || !st_facilities::Util::fileExists(pointingHistory));
   std::string sc_table = m_pars["sctable"];
   observationSim::ScDataContainer scData(prefix + "_scData", sc_table,
                                          nMaxRows, writeScData, &m_pars);
   scData.setAppName("gtobssim");
   scData.setVersion(getVersion());
   observationSim::Spacecraft * spacecraft(0);
   if (writeScData) {
      spacecraft = new observationSim::LatSc();
   } else {
      spacecraft = new observationSim::LatSc(pointingHistory);
   }
   double frac = m_pars["ltfrac"];
   spacecraft->setLivetimeFrac(frac);
   if (m_pars["nevents"]) {
      m_formatter->info() << "Generating " << m_count 
                          << " events...." << std::endl;
      m_simulator->generateEvents(static_cast<long>(m_count), events, 
                                  scData, m_respPtrs, spacecraft);
   } else {
      m_formatter->info() << "Generating events for a simulation time of "
                          << m_count << " seconds...." << std::endl;
      m_simulator->generateEvents(m_count, events, scData, m_respPtrs, 
                                  spacecraft);
   }

   if (writeScData) {
// Pad with one more row of ScData.
      double time = scData.simTime() + 30.;
      scData.addScData(time, spacecraft);
   }

   saveEventIds(events);
}

void ObsSim::
saveEventIds(const observationSim::EventContainer & events) const {
   typedef observationSim::EventContainer::SourceSummary srcSummary_t;
   typedef std::map<std::string, srcSummary_t> id_map_t;

   const id_map_t & eventIds = events.eventIds();

// sort by ID number
   unsigned int nsrcs = eventIds.size();
   std::vector<int> idnums(nsrcs);
   std::vector<std::string> ids(nsrcs);
   std::vector<unsigned long> incidents(nsrcs);
   std::vector<unsigned long> accepteds(nsrcs);
   id_map_t::const_iterator eventId = eventIds.begin();
   for (size_t idnum=0 ; eventId != eventIds.end(); ++eventId, idnum++) {
      ids.at(idnum) = eventId->first;
      idnums.at(idnum) = eventId->second.id;
      incidents.at(idnum) = eventId->second.incidentNum;
      accepteds.at(idnum) = eventId->second.acceptedNum;
   }
   
   std::string event_id_file = m_pars["evroot"];
   event_id_file += "_srcIds.txt";
   std::ofstream outputFile(event_id_file.c_str());
   for (unsigned int i = 0; i < nsrcs; i++) {
      outputFile << idnums.at(i) << "  "
                 << ids.at(i) << "  "
                 << incidents.at(i) << "  "
                 << accepteds.at(i) << "\n";
   }
   outputFile.close();
}

double ObsSim::maxEffArea() const {
   if (m_respPtrs.empty()) {
      double effArea = m_pars["area"];
      return effArea;
   }
   double total(0);
   for (size_t i=0; i < m_respPtrs.size(); i++) {
      total += m_respPtrs.at(i)->aeff()->upperLimit();
   }
   if (m_respPtrs.front()->efficiencyFactor()) {
// Provide head room for efficiency factor corrections at large
// ltfrac.
      total *= 1.5;
   }
   return total/1e4;
}

void ObsSim::get_tstart(std::string scfile, const std::string & sctable) {
   facilities::Util::expandEnvVar(&scfile);
   std::auto_ptr<const tip::Table>
      sc_data(tip::IFileSvc::instance().readTable(scfile, sctable));
/// TSTART from the Fermi astroserver is unreliable.  Use START of
/// first entry instead.
//    const tip::Header & header(sc_data->getHeader());
//    header["TSTART"].get(m_tstart);
   tip::Table::ConstIterator it = sc_data->begin();
   tip::ConstTableRecord & row = *it;
   row["start"].get(m_tstart);
}
