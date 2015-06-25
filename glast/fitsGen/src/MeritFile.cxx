/**
 * @file MeritFile.cxx
 * @brief Implementation for merit tuple file abstraction using tip.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/fitsGen/src/MeritFile.cxx,v 1.11 2010/07/21 04:43:13 jchiang Exp $
 */

#include <algorithm>

#include "tip/IFileSvc.h"

#include "fitsGen/MeritFile.h"

namespace fitsGen {

MeritFile::MeritFile(const std::string & meritfile,
                     const std::string & tree,
                     const std::string & filter)
   : m_table(tip::IFileSvc::instance().readTable(meritfile, tree, filter)),
     m_it(m_table->begin()),
     m_row(*m_it),
     m_nrows(m_table->getNumRecords()),
     m_haveTime(true) {
   const std::vector<std::string> & validFields(m_table->getValidFields());
   if (std::find(validFields.begin(), validFields.end(), "EvtElapsedTime") 
       == validFields.end()) {
      m_haveTime = false;
   } else {
      m_tstart = m_row["EvtElapsedTime"].get();
      m_it = end();
      --m_it;
      m_tstop = m_row["EvtElapsedTime"].get();
      m_it = begin();
   }
}

MeritFile::~MeritFile() {
   delete m_table;
}

void MeritFile::next() {
   ++m_it;
}

void MeritFile::prev() {
   --m_it;
}

double MeritFile::operator[](const std::string & fieldname) const {
   if (fieldname == "EvtElapsedTime" && !m_haveTime) {
      return 0;
   }
   return m_row[fieldname].get();
}

tip::Table::ConstIterator MeritFile::begin() const {
   return m_table->begin();
}

tip::Table::ConstIterator MeritFile::end() const {
   return m_table->end();
}

tip::Table::ConstIterator & MeritFile::itor() {
   return m_it;
}

void MeritFile::setStartStop(double tstart, double tstop) {
   m_tstart = tstart;
   m_tstop = tstop;
}

short int MeritFile::conversionType() const {
   if (17 - m_row["Tkr1FirstLayer"].get() < 11.5) { // Front converting
      return 0;
   } else { // Back converting
      return 1;
   }
}

} // namespace fitsGen
