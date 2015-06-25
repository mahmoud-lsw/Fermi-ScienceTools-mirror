/**
 * @file Cuts.cxx
 * @brief Handle data selections and DSS keyword management.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/src/Cuts.cxx,v 1.68.2.1 2015/04/23 20:50:49 jchiang Exp $
 */

#include <cctype>
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "irfUtil/Util.h"

#include "dataSubselector/BitMaskCut.h"
#include "dataSubselector/Cuts.h"
#include "dataSubselector/GtiCut.h"
#include "dataSubselector/RangeCut.h"
#include "dataSubselector/SkyConeCut.h"
#include "dataSubselector/VersionCut.h"

namespace {
void toUpper(std::string & name) {
   for (std::string::iterator it = name.begin(); it != name.end(); ++it) {
      *it = std::toupper(*it);
   }
}

dataSubselector::RangeCut * 
mergeRangeCuts(const std::vector<dataSubselector::RangeCut *> & cuts) {
   double minValue;
   double maxValue;
   bool haveMin(false);
   bool haveMax(false);
   for (size_t i = 0; i < cuts.size(); i++) {
      if (cuts.at(i)->intervalType() == dataSubselector::RangeCut::CLOSED ||
          cuts.at(i)->intervalType() == dataSubselector::RangeCut::MINONLY) {
         if (!haveMin) {
            minValue = cuts.at(i)->minVal();
            haveMin = true;
         } else if (cuts.at(i)->minVal() > minValue) {
            minValue = cuts.at(i)->minVal();
         }
      }
      if (cuts.at(i)->intervalType() == dataSubselector::RangeCut::CLOSED ||
          cuts.at(i)->intervalType() == dataSubselector::RangeCut::MAXONLY) {
         if (!haveMax) {
            maxValue = cuts.at(i)->maxVal();
            haveMax = true;
         } else if (cuts.at(i)->maxVal() < maxValue) {
            maxValue = cuts.at(i)->maxVal();
         }
      }
   }
   dataSubselector::RangeCut::IntervalType 
      type(dataSubselector::RangeCut::CLOSED);
   if (!haveMin) {
      type = dataSubselector::RangeCut::MAXONLY;
   } else if (!haveMax) {
      type = dataSubselector::RangeCut::MINONLY;
   }
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(cuts.at(0)->colname(), "[]", tokens);
   return new dataSubselector::RangeCut(tokens.at(0),
                                        cuts.at(0)->unit(),
                                        minValue, maxValue, type,
                                        cuts.at(0)->index());
}

void joinPaths(const std::string & components,
               std::string & final_path) {
   std::vector<std::string> tokens; 
   facilities::Util::stringTokenize(components, " ", tokens);
   if (tokens.size() == 0) {
      throw std::runtime_error("dataSubselector::Cuts: joinPaths: "
                               "empty string input");
   }
   if (tokens.size() == 1) {
      final_path = tokens[0];
      return;
   }
   final_path = facilities::commonUtilities::joinPath(tokens[0], tokens[1]);
   for (size_t i(2); i < tokens.size(); i++) {
      final_path = facilities::commonUtilities::joinPath(final_path,tokens[i]);
   }
   return;
}

} // anonymous namespace

namespace dataSubselector {

Cuts::Cuts(const std::vector<std::string> & eventFiles,
           const std::string & extname, bool check_columns,
           bool skipTimeRangeCuts, bool skipEventClassCuts) 
   : m_irfName("NONE"), m_post_P7(false) {
   std::vector<Cuts> my_cuts;
   for (size_t i = 0; i < eventFiles.size(); i++) {
      my_cuts.push_back(Cuts(eventFiles.at(i), extname, check_columns,
                             skipTimeRangeCuts, skipEventClassCuts));
      if (i > 0) {
         if (!my_cuts.front().compareWithoutGtis(my_cuts.back())) {
            std::ostringstream message;
            message << "DSS keywords in " << eventFiles.at(i)
                    << " do not match those in " << eventFiles.front();
            throw std::runtime_error(message.str());
         }
      }
   }
   if (!my_cuts.empty()) {
      *this = mergeGtis(my_cuts);
   }
   set_irfName(eventFiles.at(0), extname);
}

Cuts Cuts::mergeGtis(std::vector<Cuts> & cuts_vector) {
// Gather non-Gti cuts.
   dataSubselector::Cuts my_cuts;
   const dataSubselector::Cuts & firstCuts(cuts_vector.front());
   for (size_t i = 0; i < firstCuts.size(); i++) {
      if (firstCuts[i].type() != "GTI") {
         my_cuts.addCut(firstCuts[i]);
      }
   }
   
// Merge all of the GTIs into one, taking the union of the intervals.
   dataSubselector::Gti merged_gti;
   std::vector<const dataSubselector::GtiCut *> gtiCuts;
   for (size_t i = 0; i < cuts_vector.size(); i++) {
      cuts_vector.at(i).getGtiCuts(gtiCuts);
      for (size_t j = 0; j < gtiCuts.size(); j++) {
         if (i == 0 && j == 0) {
            merged_gti = gtiCuts.at(j)->gti();
         } else {
            merged_gti = merged_gti | gtiCuts.at(j)->gti();
         }
      }
   }

   if (merged_gti.getNumIntervals() > 0) {
      my_cuts.addGtiCut(merged_gti);
   }
   return my_cuts;
}

void Cuts::getGtiCuts(std::vector<const GtiCut *> & gtiCuts) {
   gtiCuts.clear();
   for (size_t i = 0; i < m_cuts.size(); i++) {
      if (m_cuts.at(i)->type() == "GTI") {
         gtiCuts.push_back(dynamic_cast<GtiCut *>(m_cuts.at(i)));
      }
   }
}

Cuts::Cuts(const std::string & eventFile, const std::string & extname,
           bool check_columns, bool skipTimeRangeCuts,
           bool skipEventClassCuts) 
   : m_irfName("NONE") {
   /// Read in validity masks for Pass 8 event type and event class
   /// selections.
   if (BitMaskCut::evclassValidityMasks() == 0 ||
       BitMaskCut::evtypeValidityMasks() == 0) {
      std::string evclassPath;
      ::joinPaths("data glast lat bcf valid_evclass_selections.txt",
                  evclassPath);
      std::string evtypePath;
      ::joinPaths("data glast lat bcf valid_evtype_selections.txt",
                  evtypePath);
      evclassPath = facilities::commonUtilities::joinPath(
         st_facilities::Environment::getEnv("CALDB"), evclassPath);
      evtypePath = facilities::commonUtilities::joinPath(
         st_facilities::Environment::getEnv("CALDB"), evtypePath);
      BitMaskCut::setValidityMasks(evclassPath, evtypePath);
   }

   const tip::Extension * ext(0);
   try {
      ext = tip::IFileSvc::instance().readTable(eventFile, extname);
   } catch (tip::TipException & eObj) {
      if (!st_facilities::Util::expectedException(eObj, "HDU is not a")) {
         throw;
      }
      ext = tip::IFileSvc::instance().readImage(eventFile, extname);
   }

   std::vector<std::string> colnames;
   if (check_columns) {
      const tip::Table * table 
         = dynamic_cast<tip::Table *>(const_cast<tip::Extension *>(ext));
      colnames = table->getValidFields();
// FITS column names are in CAPS, not lowercase, so undo what tip has wrought
      for (size_t i = 0; i < colnames.size(); i++) {
         ::toUpper(colnames[i]);
      }
   }

   const tip::Header & header = ext->getHeader();

// NB: The .get(...) method does not work for unsigned int arguments.
   int nkeys;
   header["NDSKEYS"].get(nkeys);

   std::string type, unit, value, ref("");

   for (int keynum = 1; keynum <= nkeys; keynum++) {
      std::ostringstream key1, key2, key3, key4;
      key1 << "DSTYP" << keynum;
      header[key1.str()].get(type);
      key2 << "DSUNI" << keynum;
      header[key2.str()].get(unit);
      key3 << "DSVAL" << keynum;
      header[key3.str()].get(value);
      ::toUpper(value);

      if (value == "TABLE") {
         key4 << "DSREF" << keynum;
         header[key4.str()].get(ref);
      }
      std::string colname;
      unsigned int indx = parseColname(type, colname);
      if (value.find("CIRCLE") == 0) {
         m_cuts.push_back(new SkyConeCut(type, unit, value));
      } else if (type == "TIME" && value == "TABLE") {
         m_cuts.push_back(new GtiCut(eventFile));
      } else if (type.substr(0, 8) == "BIT_MASK") {
         std::vector<std::string> tokens;
         facilities::Util::stringTokenize(type, "(),", tokens);
         unsigned int mask = std::atoi(tokens[2].c_str());
         if (tokens.size() == 4) {
            // The third position in the BIT_MASK arg list is pass_ver.
            if (!BitMaskCut::post_P7(tokens[3])) {
               // For pre-Pass 8 data, the value of mask is the bit
               // position so do the bit shift to generate the mask.
               mask = 1 << mask;
            }
            m_cuts.push_back(new BitMaskCut(tokens[1], mask, tokens[3]));
         } else {
            // This is also (only) pre-Pass 8 and probably cannot
            // occur anymore.
            mask = 1 << mask;
            m_cuts.push_back(new BitMaskCut(tokens[1], mask));
         }
      } else if (type.length() >= 7 &&
                 type.substr(type.length()-7, 7) == "VERSION") {
         m_cuts.push_back(new VersionCut(colname, value));
      } else if ( (!check_columns || 
                   std::find(colnames.begin(), colnames.end(), colname) 
                   != colnames.end())
                  && value != "TABLE" ) {
         if ((type != "TIME" || !skipTimeRangeCuts) &&
             (type != "EVENT_CLASS" || !skipEventClassCuts)) {
            m_cuts.push_back(new RangeCut(colname, unit, value, indx));
         }
      } else {
         std::ostringstream message;
         message << "FITS extension contains an unrecognized DSS filtering "
                 << "scheme.\n"
                 << key1.str() << " = " << type << "\n"
                 << key2.str() << " = " << unit << "\n"
                 << key3.str() << " = " << value << "\n";
         if (value == "TABLE") {
            message << key4.str() << " = " << ref << "\n";
         }
         throw std::runtime_error(message.str());
      }
   }
   delete ext;
   set_irfName(eventFile, extname);
   read_pass_ver(eventFile, extname);
}

unsigned int Cuts::parseColname(const std::string & colname,
                                std::string & col) const {
   if (colname.find("[") == std::string::npos) {
      col = colname;
      return 0;
   }
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(colname, "[]", tokens);
   col = tokens.at(0);
   return std::atoi(tokens.at(1).c_str());
}

unsigned int Cuts::find(const CutBase * cut) const {
   for (unsigned int i = 0; i < size(); i++) {
      if (*cut == *m_cuts[i]) {
         return i;
      }
   }
   return size();
}

bool Cuts::hasCut(const CutBase * newCut) const {
   for (unsigned int i = 0; i < size(); i++) {
      if (*newCut == *m_cuts[i]) {
         return true;
      }
   }
   return false;
}

Cuts::Cuts(const Cuts & rhs) 
   : m_irfName(rhs.m_irfName),
     m_pass_ver(rhs.m_pass_ver) {
   m_cuts.reserve(rhs.size());
   for (unsigned int i = 0; i < rhs.size(); i++) {
      m_cuts.push_back(rhs.m_cuts[i]->clone());
   }
}

Cuts::~Cuts() {
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      delete m_cuts[i];
   }
}

Cuts & Cuts::operator=(const Cuts & rhs) {
   if (*this != rhs) {
      std::vector<CutBase *>::reverse_iterator cut = m_cuts.rbegin();
      for ( ; cut != m_cuts.rend(); cut++) {
         delete *cut;
      }
      m_cuts.clear();
      for (unsigned int i = 0; i < rhs.size(); i++) {
         m_cuts.push_back(rhs.m_cuts.at(i)->clone());
      }
   }
   return *this;
}

bool Cuts::accept(tip::ConstTableRecord & row) const {
   bool ok(true);
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      ok = ok && m_cuts[i]->accept(row);
   }
   return ok;
}

bool Cuts::accept(const std::map<std::string, double> & params) const {
   bool ok(true);
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      ok = ok && m_cuts[i]->accept(params);
   }
   return ok;
}

unsigned int Cuts::addRangeCut(const std::string & colname,
                               const std::string & unit,
                               double minVal, double maxVal, 
                               RangeCut::IntervalType type,
                               unsigned int indx) {
   return addCut(new RangeCut(colname, unit, minVal, maxVal, type, indx));
}

unsigned int Cuts::addGtiCut(const tip::Table & table) {
   return addCut(new GtiCut(table));
}

unsigned int Cuts::addGtiCut(const Gti & gti) {
   return addCut(new GtiCut(gti));
}

unsigned int Cuts::addSkyConeCut(double ra, double dec, double radius) {
   return addCut(new SkyConeCut(ra, dec, radius));
}

unsigned int Cuts::addBitMaskCut(const std::string & colname,
                                 unsigned int mask,
                                 const std::string & pass_ver) {
   return addCut(new BitMaskCut(colname, mask, pass_ver));
}

unsigned int Cuts::addVersionCut(const std::string & colname,
                                 const std::string & version) {
   return addCut(new VersionCut(colname, version));
}

unsigned int Cuts::addCut(CutBase * newCut) {
   if (hasCut(newCut)) {
      delete newCut;
   } else {
      for (unsigned int j = 0; j != size(); j++) {
         if (newCut->supercedes(*(m_cuts[j]))) {
            delete m_cuts[j];
            m_cuts[j] = newCut;
            return size();
         }
         if (m_cuts[j]->supercedes(*newCut)) {
            delete newCut;
            return size();
         }
      }
      m_cuts.push_back(newCut);
   }
   return size();
}

unsigned int Cuts::mergeRangeCuts() {
   std::vector<RangeCut *> rangeCuts;
   for (size_t j = 0; j < m_cuts.size(); j++) {
      if (m_cuts.at(j)->type() == "range") {
         RangeCut * rangeCut(dynamic_cast<RangeCut *>(m_cuts.at(j)));
         rangeCuts.push_back(rangeCut);
      }
   }

   std::map<std::string, int> colnames;
   for (size_t j = 0; j < rangeCuts.size(); j++) {
      colnames[rangeCuts.at(j)->colname()] = 1;
   }

   for (std::map<std::string, int>::const_iterator colname(colnames.begin());
        colname != colnames.end(); ++colname) {
      removeRangeCuts(colname->first, rangeCuts);
      m_cuts.push_back(::mergeRangeCuts(rangeCuts));
      for (size_t i = 0; i < rangeCuts.size(); i++) {
         delete rangeCuts.at(i);
      }
   }

   return m_cuts.size();
}

unsigned int Cuts::removeVersionCut(const std::string & colname) {
   std::vector<CutBase *> held_cuts;
   for (size_t i(0); i < m_cuts.size(); i++) {
      if (m_cuts.at(i)->type() == "version") {
         VersionCut * versionCut(dynamic_cast<VersionCut *>(m_cuts.at(i)));
         if (versionCut->colname() == colname) {
            delete m_cuts.at(i);
            continue;
         }
      }
      held_cuts.push_back(m_cuts.at(i));
   }
   m_cuts = held_cuts;
   return m_cuts.size();
}

unsigned int Cuts::removeRangeCuts(const std::string & colname,
                                   std::vector<RangeCut *> & removedCuts) {
   removedCuts.clear();
   for (size_t j = 0; j < m_cuts.size(); j++) {
      if (m_cuts.at(j)->type() == "range") {
         RangeCut * rangeCut(dynamic_cast<RangeCut *>(m_cuts.at(j)));
         if (rangeCut->colname() == colname) {
            removedCuts.push_back(rangeCut);
            m_cuts.at(j) = 0;
         }
      }
   }
   std::vector<CutBase *> held_cuts;
   for (size_t j = 0; j < m_cuts.size(); j++) {
      if (m_cuts.at(j)) {
         held_cuts.push_back(m_cuts.at(j));
      }
   }
   m_cuts = held_cuts;
   return m_cuts.size();
}

void Cuts::writeDssKeywords(tip::Header & header) const {
   removeDssKeywords(header);
   int ndskeys(0);
   for (size_t i = 0; i < m_cuts.size(); i++) {
      if (!isTimeCut(*m_cuts.at(i)) || m_cuts.at(i)->type() == "GTI") {
         ndskeys++;
         m_cuts[i]->writeDssKeywords(header, ndskeys);
      }
   }
   header["NDSKEYS"].set(ndskeys);
}

void Cuts::writeDssTimeKeywords(tip::Header & header) const {
   removeDssKeywords(header);

   std::vector<CutBase *> my_time_cuts;
   for (size_t i = 0; i < m_cuts.size(); i++) {
      if (isTimeCut(*m_cuts.at(i))) {
         my_time_cuts.push_back(m_cuts.at(i));
      }
   }

   int ndskeys = my_time_cuts.size();
   header["NDSKEYS"].set(ndskeys);
   for (unsigned int i = 0; i < my_time_cuts.size(); i++) {
      my_time_cuts.at(i)->writeDssKeywords(header, i + 1);
   }
}

bool Cuts::isTimeCut(const CutBase & cut) {
   if (cut.type() == "GTI" || 
       (cut.type() == "range" && 
        dynamic_cast<RangeCut &>(const_cast<CutBase &>(cut)).colname() 
        == "TIME")) {
      return true;
   }
   return false;
}

void Cuts::checkIrfs(const std::string & infile, 
                     const std::string & extname,
                     const std::string & irfs) {
   bool check_columns;
   Cuts my_cuts(infile, extname, check_columns=false);
   if (my_cuts.irfName() != "NONE" && my_cuts.irfName() != irfs) {
      st_stream::StreamFormatter formatter("dataSubselector::Cuts",
                                           "checkIrfs", 2);
      formatter.warn() << "\nWARNING:\n"
                       << "IRF version mismatch detected. "
                       << "IRF version in HEADER: " 
                       << my_cuts.irfName() << ", "
                       << "IRF version provided (by CALDB/on command line): "
                       << irfs <<  std::endl;
   }
}

void Cuts::removeDssKeywords(tip::Header & header) const {
   int ndskeys(0);
   try {
      header["NDSKEYS"].get(ndskeys);
      char * dskeys[] = {const_cast<char *>("DSTYP"), 
                         const_cast<char *>("DSUNI"), 
                         const_cast<char *>("DSVAL"),
                         const_cast<char *>("DSREF")};
      for (int i = 0; i < ndskeys; i++) {
         for (int j = 0; j < 4; j++) {
            std::ostringstream keyname;
            keyname << dskeys[j] << i+1;
            tip::Header::Iterator keyword = header.find(keyname.str());
            if (keyword != header.end()) {
               header.erase(keyword);
            }
         }
      }
      header.erase("NDSKEYS");
   } catch (tip::TipException & eObj) {
      if (!st_facilities::Util::expectedException(eObj, "keyword not found")) {
         throw;
      }
   }
}

void Cuts::writeGtiExtension(const std::string & filename) const {
   const dataSubselector::Gti * gti(0);
   for (unsigned int i = 0; i < size(); i++) {
      if (m_cuts[i]->type() == "GTI") {
         gti = &(dynamic_cast<GtiCut *>(m_cuts[i])->gti());
         break;  // assume there is at most one GTI
      }
   }
   if (gti) {
      gti->writeExtension(filename);
   }
}

bool Cuts::operator==(const Cuts & rhs) const {
   if (size() != rhs.size()) {
      return false;
   }
   for (unsigned int i = 0; i < size(); i++) {
      unsigned int place = find(rhs.m_cuts.at(i));
      if (place == size()){
         return false;
      }
   }
   return true;
}

bool Cuts::compareWithoutGtis(const Cuts & rhs) const {
   if (size() != rhs.size()) {
      return false;
   }
   for (unsigned int i = 0; i < size(); i++) {
      if (rhs.m_cuts.at(i)->type() != "GTI") {
         unsigned int place = find(rhs.m_cuts.at(i));
         if (place == size()){
            return false;
         }
      }
   }
   return true;
}

void Cuts::writeCuts(std::ostream & stream, bool suppressGtis) const {
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      if (!suppressGtis || m_cuts.at(i)->type() != "GTI") {
         m_cuts.at(i)->writeCut(stream, i + 1);
      } else {
         stream << "DSTYP" << i+1 << ": TIME\n"
                << "DSUNI" << i+1 << ": s\n"
                << "DSVAL" << i+1 << ": TABLE\n"
                << "DSREF" << i+1 << ": :GTI\n\n"
                << "GTIs: (suppressed)\n\n";
      }
   }
}

std::string Cuts::filterString() const {
   std::string filter("");
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      std::string my_filter = m_cuts.at(i)->filterString();
      if (my_filter != "") {
         if (filter != "") {
            filter += " && ";
         }
         filter += my_filter;
      }
   }
   return filter;
}

const std::string & Cuts::irfName() const {
   return m_irfName;
}

std::string Cuts::CALDB_implied_irfs() const {
   std::map<std::string, unsigned int> irfs;
   read_bitmask_mapping(irfs);
   // Test against names in the EVENT_CLASS column of the
   // BITMASK_MAPPING extension of irf_index.fits.  
   //
   // Strip partition name before testing.
   std::string test_irfName(m_irfName);
   std::string::size_type pos;
   if ((pos = test_irfName.find(" (")) != std::string::npos) {
      test_irfName = test_irfName.substr(0, pos);
   }
   if (test_irfName != "NONE" && irfs.find(test_irfName) == irfs.end()) {
      throw std::runtime_error("Invalid IRF name: " + test_irfName);
   }
   const BitMaskCut * my_bitmask_cut(bitMaskCut("EVENT_CLASS"));
   if (my_bitmask_cut == 0) {
      throw std::runtime_error("No EVENT_CLASS bitmask cut in input file, so "
                               "cannot infer most recent IRFs from CALDB.");
   }
   unsigned int mask(my_bitmask_cut->mask());
   delete my_bitmask_cut;
   std::map<std::string, unsigned int>::const_iterator it(irfs.begin());
   std::string irfs_name("");
   unsigned int irf_ver_num(0);
   for ( ; it != irfs.end(); ++it) {
      if (it->second != mask) {
         continue;
      }
      std::string pass_ver;
      std::string irf_ver;
      extract_irf_versions(it->first, pass_ver, irf_ver);
      unsigned int candidate_irf_ver_num(std::atoi(irf_ver.substr(1).c_str()));
      if (pass_ver == m_pass_ver && 
          (irfs_name == "" || candidate_irf_ver_num > irf_ver_num)) {
         irfs_name = it->first;
         irf_ver_num = candidate_irf_ver_num;
      }
   }
   append_event_type_partition(irfs_name);
   return irfs_name;
}

void Cuts::append_event_type_partition(std::string & irfs_name) const {
   /// Infer event_type partition.
   const BitMaskCut * event_type_cut(bitMaskCut("EVENT_TYPE"));
   if (event_type_cut) {
      typedef std::map<std::string, std::pair<unsigned int, std::string> > 
         EventTypeMapping_t;
      EventTypeMapping_t evtype_mapping;
      irfUtil::Util::get_event_type_mapping(irfs_name, evtype_mapping);

      unsigned int bit(static_cast<int>(std::log(event_type_cut->mask())/
                                        std::log(2)));
      for (EventTypeMapping_t::const_iterator it(evtype_mapping.begin());
           it != evtype_mapping.end(); ++it) {
         if (bit == it->second.first && it->second.second != "none") {
            irfs_name += (" (" + it->second.second + ")");
            break;
         }
      }
   }
}

void Cuts::
read_bitmask_mapping(std::map<std::string, unsigned int> & irfs) {
   irfs.clear();
   std::string sub_path;
   ::joinPaths("data glast lat bcf irf_index.fits", sub_path);
   std::string irf_index = facilities::commonUtilities::joinPath(
      st_facilities::Environment::getEnv("CALDB"), sub_path);
   const tip::Table * irf_map 
      = tip::IFileSvc::instance().readTable(irf_index, "BITMASK_MAPPING");
   tip::Table::ConstIterator it(irf_map->begin());
   tip::ConstTableRecord & row = *it;
   for ( ; it != irf_map->end(); ++it) {
      std::string event_class;
      int bitpos;
      row["event_class"].get(event_class);
      row["bitposition"].get(bitpos);
      irfs[event_class] = 1 << bitpos;
   }
   delete irf_map;
}

void Cuts::read_pass_ver(const std::string & infile, 
                         const std::string & ext) {
// Set default value.
   m_pass_ver == "NONE";

// Attempt to read m_pass_ver from PASS_VER keyword.
/// @todo Improve error handling when PASS_VER does not exist or does not
/// have "P#V#" format.
   std::auto_ptr<const tip::Extension> 
      events(tip::IFileSvc::instance().readExtension(infile, ext));
   try {
      events->getHeader()["PASS_VER"].get(m_pass_ver);
   } catch (tip::TipException & eObj) {
      if (st_facilities::Util::expectedException(eObj,"Cannot read keyword")) {
         // Look for pass_ver from BIT_MASK cut.
         if (bitMaskCut()) {
            m_pass_ver = bitMaskCut()->pass_ver();
            return;
         }
      }
   }
   return;
}

void Cuts::set_irfName(const std::string & infile, 
                       const std::string & ext) {
   for (size_t i(0); i < m_cuts.size(); i++) {
      VersionCut * version_cut 
         = dynamic_cast<VersionCut *>(const_cast<CutBase *>(m_cuts[i]));
      if (version_cut && version_cut->colname() == "IRF_VERSION") {
         m_irfName = version_cut->version();
//         continue;
         break;
      }
/// @todo Need to determine if there is any context where adding the
/// FRONT/BACK qualifier is needed, since CONVTYPE cuts are included
/// in DSS keywords already and read in by tools like gtexpmap.
//
//       RangeCut * convtype_cut 
//          = dynamic_cast<RangeCut *>(const_cast<CutBase *>(m_cuts[i]));
//       if (convtype_cut && convtype_cut->colname() == "CONVERSION_TYPE") {
//          if (convtype_cut->minVal() == 0 && 
//              convtype_cut->maxVal() == 0) {
//             m_irfName += "::FRONT";
//          } else if (convtype_cut->minVal() == 1 && 
//                     convtype_cut->maxVal() == 1) {
//             m_irfName += "::BACK";
//          }
//       }
   }
}

void Cuts::extract_irf_versions(const std::string & irf_name,
                                std::string & pass_ver,
                                std::string & irf_ver) {
   size_t evlen(irf_name.length());
   size_t v_pos(irf_name.find_last_of('V'));
   if (v_pos == std::string::npos) {
      pass_ver = "NONE";
      irf_ver = "NONE";
      return;
   } else if (irf_name.substr(0, 2) == "P7" &&
              irf_name.substr(v_pos, evlen - v_pos) == "V6") {
      // Handle Pass 7, V6 irfs as a special case.
      pass_ver = "P7V6";
      irf_ver = "V6";
      return;
   }
   // The PASS_VER value in the irf name should be set off from 
   // the remainder of the irf name by a "_" delimiter.
   size_t pv_delim(irf_name.find_first_of('_'));
   pass_ver = irf_name.substr(0, pv_delim);
   irf_ver = irf_name.substr(v_pos, evlen - v_pos);
}

void Cuts::setIrfs(const std::string & irfName) {
   if (bitMaskCut()) {
      throw std::runtime_error("dataSubselector::Cuts::setIrfs: "
                               "Bit mask cut already set");
   }
   m_irfName = irfName;
   size_t delim_pos = irfName.find_last_of(":");
   if (delim_pos != std::string::npos) { 
      // Have a CONVERSION_TYPE selection, so check against existing
      // selection.
      std::string section(irfName.substr(delim_pos - 1,
                                         irfName.length() - delim_pos + 1));
      RangeCut * convtype_cut(0);
      for (size_t i(0); i < m_cuts.size(); i++) {
         convtype_cut 
            = dynamic_cast<RangeCut *>(const_cast<CutBase *>(m_cuts[i]));
         if (convtype_cut && convtype_cut->colname() == "CONVERSION_TYPE") {
            if ((convtype_cut->minVal() == 0 && 
                 convtype_cut->maxVal() == 0 &&
                 section != "::FRONT") ||
                (convtype_cut->minVal() == 1 && 
                 convtype_cut->maxVal() == 1 &&
                 section != "::BACK")) {
               throw std::runtime_error("dataSubselector::Cuts::setIrfs: "
                                        "Inconsistent FRONT/BACK selection");
            }
         }
      }
      if (!convtype_cut) {
         // Add a CONVERSION_TYPE selection based in irfName.
         if (section == "::FRONT") {
            addRangeCut("CONVERSION_TYPE", "dimensionless", 0, 0);
         } else if (section == "::BACK") {
            addRangeCut("CONVERSION_TYPE", "dimensionless", 1, 1);
         } else {
            throw std::runtime_error("invalid section: " + section);
         }
      }
   }

   std::string event_class(irfName.substr(0, delim_pos - 1));
   std::string irf_ver;
   extract_irf_versions(event_class, m_pass_ver, irf_ver);
   m_post_P7 = BitMaskCut::post_P7(m_pass_ver);
   
   std::map<std::string, unsigned int> irfs;
   try {
      read_bitmask_mapping(irfs);
      std::map<std::string, unsigned int>::const_iterator 
         it = irfs.find(event_class);
      if (it != irfs.end()) {
         unsigned int mask(it->second);
         addBitMaskCut("EVENT_CLASS", mask, m_pass_ver);
      }
   } catch (std::runtime_error & eObj) {
      if (!st_facilities::Util::expectedException(eObj,
                                                  "read_bitmask_mapping")) {
         throw;
      }
   }
   addVersionCut("IRF_VERSION", irfName);
}

BitMaskCut * Cuts::bitMaskCut(const std::string & colname) const {
   for (size_t i(0); i < m_cuts.size(); i++) {
      const BitMaskCut * bit_mask_cut 
         = dynamic_cast<const BitMaskCut *>(m_cuts[i]);
      if (bit_mask_cut && bit_mask_cut->colname() == colname) {
         return new BitMaskCut(*bit_mask_cut);
      }
   }
   return 0;
}

std::vector<BitMaskCut *> Cuts::bitMaskCuts() const {
   std::vector<BitMaskCut *> my_bitMaskCuts;
   for (size_t i(0); i < m_cuts.size(); i++) {
      const BitMaskCut * bit_mask_cut
         = dynamic_cast<const BitMaskCut *>(m_cuts[i]);
      if (bit_mask_cut) {
         my_bitMaskCuts.push_back(new BitMaskCut(*bit_mask_cut));
      }
   }
   return my_bitMaskCuts;
}

RangeCut * Cuts::conversionTypeCut() const {
   for (size_t i(0); i < m_cuts.size(); i++) {
      RangeCut * range_cut 
         = dynamic_cast<RangeCut *>(const_cast<CutBase *>(m_cuts[i]));
      if (range_cut && range_cut->colname() == "CONVERSION_TYPE") {
         return range_cut;
      }
   }
   return 0;
}

void Cuts::setBitMaskCut(BitMaskCut * candidateCut) {
   if (!candidateCut) {
      // Do nothing with a null pointer.
      return;
   }
   // Delete any existing BitMaskCuts operating on the same column.
   std::vector<CutBase *> my_cuts;
   for (size_t i(0); i < m_cuts.size(); i++) {
      BitMaskCut * currentCut(dynamic_cast<BitMaskCut *>(m_cuts[i]));
      if (currentCut && (currentCut->colname() == candidateCut->colname())) {
         delete m_cuts[i];
      } else {
         my_cuts.push_back(m_cuts[i]);
      }
   }
   // Add the candidate cut.
   my_cuts.push_back(candidateCut);
   m_cuts = my_cuts;
}

} // namespace dataSubselector
