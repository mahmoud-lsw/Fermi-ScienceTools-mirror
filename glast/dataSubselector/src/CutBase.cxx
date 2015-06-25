/**
 * @file CutBase.cxx
 * @brief Implementation for CutBase methods.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/src/CutBase.cxx,v 1.3 2005/06/01 19:49:33 jchiang Exp $
 */

#include <sstream>

#include "tip/Header.h"

#include "dataSubselector/CutBase.h"

namespace dataSubselector {

bool CutBase::operator==(const CutBase & rhs) const {
   return m_type == rhs.m_type && this->equals(rhs);
}

void CutBase::writeCut(std::ostream & stream, unsigned int keynum) const {
   std::string type, unit, value, ref("");
   getKeyValues(type, unit, value, ref);
   stream << "DSTYP" << keynum << ": " << type << "\n"
          << "DSUNI" << keynum << ": " << unit << "\n"
          << "DSVAL" << keynum << ": " << value << "\n";
   if (ref != "") {
      stream << "DSREF" << keynum << ": " << ref << "\n";
   }
   stream << std::endl;
}

void CutBase::writeDssKeywords(tip::Header & header,
                               unsigned int keynum) const {
   std::string type, unit, value, ref("");
   getKeyValues(type, unit, value, ref);
   writeDssKeywords(header, keynum, type, unit, value, ref);
}

void CutBase::writeDssKeywords(tip::Header & header, 
                               unsigned int keynum,
                               const std::string & type, 
                               const std::string & unit,
                               const std::string & value,
                               const std::string & ref) const {
   std::ostringstream key1, key2, key3;
   key1 << "DSTYP" << keynum;
   header[key1.str()].set(type);
   key2 << "DSUNI" << keynum;
   header[key2.str()].set(unit);
   key3 << "DSVAL" << keynum;
   header[key3.str()].set(value);
   if (ref != "") {
      std::ostringstream key4;
      key4 << "DSREF" << keynum;
      header[key4.str()].set(ref);
   }
}

} // namespace dataSubselector
