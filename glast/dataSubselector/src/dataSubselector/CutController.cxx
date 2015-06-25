/**
 * @file CutController.cxx
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/src/dataSubselector/CutController.cxx,v 1.32 2015/02/22 01:11:03 jchiang Exp $
 */

#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"
#include "st_app/AppParGroup.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "astro/SkyDir.h"

#include "dataSubselector/BitMaskCut.h"
#include "dataSubselector/CutController.h"
#include "dataSubselector/Gti.h"

namespace dataSubselector {

CutController * CutController::s_instance(0);

CutController * 
CutController::instance(st_app::AppParGroup & pars, 
                        const std::vector<std::string> & eventFiles,
                        const std::string & evtable) {
   if (!s_instance) {
      s_instance = new CutController(pars, eventFiles, evtable);
   }
   return s_instance;
}

void CutController::delete_instance() {
   delete s_instance;
   s_instance = 0;
}

CutController::CutController(st_app::AppParGroup & pars, 
                             const std::vector<std::string> & eventFiles,
                             const std::string & evtable) 
   : m_pars(pars), m_cuts(eventFiles, evtable, true, true), 
     m_passVer(""), m_evclsFilter("") {
   checkPassVersion(eventFiles);
   double ra = pars["ra"];
   double dec = pars["dec"];
   double radius = pars["rad"];
   double max_rad = 180.;
   if (radius < max_rad) {
      m_cuts.addSkyConeCut(ra, dec, radius);
   }
   addRangeCut("TIME", "s", pars["tmin"], pars["tmax"]);
   addRangeCut("ENERGY", "MeV", pars["emin"], pars["emax"]);
   int convtype = pars["convtype"];
   if (convtype >= 0) {
      addRangeCut("CONVERSION_TYPE", "dimensionless", convtype, convtype, 
                  0, true);
   }
   if (m_passVer == "NONE") {
      int evclsmin = pars["evclsmin"];
      int evclsmax = pars["evclsmax"];
      if (evclsmin !=0 || evclsmax != 10) {
         addRangeCut("EVENT_CLASS", "dimensionless", evclsmin, evclsmax);
      }
   } else {
      try {
         int evclass = pars["evclass"];
         if (!BitMaskCut::post_P7(m_passVer)) {
            // Pass 7 and earlier expects evclass to be the bit
            // position, so do the bit shift.
            evclass = 1 << evclass;
         }
         m_cuts.addBitMaskCut("EVENT_CLASS", evclass, m_passVer);
      } catch (const hoops::Hexception &) {
         // Assume INDEF is given as the parameter value for evclass,
         // so use default of applying no EVENT_CLASS cut.
      }
      if (BitMaskCut::post_P7(m_passVer)) {
         try {
            /// Handle any EVENT_TYPE cut supplied by the user.
            int evtype = pars["evtype"];
            m_cuts.addBitMaskCut("EVENT_TYPE", evtype, m_passVer);
         } catch (const hoops::Hexception &) {
            // Assume INDEF is given as the parameter value for evtype,
            // so use default of applying no EVENT_TYPE cut.
         }
      }
   }
   double zmin = pars["zmin"];
   double zmax = pars["zmax"];
   if (zmin > 0 || zmax < 180.) {
      addRangeCut("ZENITH_ANGLE", "deg", zmin, zmax);
   }
   double phasemin = pars["phasemin"];
   double phasemax = pars["phasemax"];
   if (phasemin !=0 || phasemax !=1) {
      addRangeCut("PULSE_PHASE", "dimensionless", phasemin, phasemax);
   }
   m_cuts.mergeRangeCuts();
}

void CutController::
checkPassVersion(const std::vector<std::string> & evfiles) {
   for (size_t i(0); i < evfiles.size(); i++) {
      const tip::Table * table = 
         tip::IFileSvc::instance().readTable(evfiles.at(i), "EVENTS");
      const tip::Header & header(table->getHeader());
      std::string passVer("NONE");
      try {
         header["PASS_VER"].get(passVer);
      } catch (tip::TipException &) {
      }
      if (m_passVer == "") {
         m_passVer = passVer;
      } else if (passVer != m_passVer) {
         delete table;
         throw std::runtime_error("PASS_VER is not consistent across "
                                  "input FT1 files");
      }
      delete table;
   }
}

bool CutController::accept(tip::ConstTableRecord & row) const {
   return m_cuts.accept(row);
}

void CutController::addRangeCut(const std::string & colname,
                                const std::string & unit,
                                double minVal, double maxVal, 
                                unsigned int indx, bool force) {
   RangeCut::IntervalType type(RangeCut::CLOSED);
   if (!force && minVal == 0 && maxVal == 0) {
      /// don't apply any range cut
      return;
   }
   if (minVal > maxVal) {
      std::ostringstream message;
      message << "minimum requested value, " << minVal 
              << ", is greater than the maximum requested, "
              << maxVal << ", for field " << colname << "\n";
      throw std::runtime_error(message.str());
   }
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(colname, "[]", tokens);
   m_cuts.addRangeCut(tokens.at(0), unit, minVal, maxVal, type, indx);
}

void CutController::updateGti(const std::string & eventFile) const {
   Gti gti(eventFile);
   for (unsigned int i = 0; i < m_cuts.size(); i++) {
      if (m_cuts[i].type() == "range") {
         const RangeCut & my_cut = 
            dynamic_cast<RangeCut &>(const_cast<CutBase &>(m_cuts[i]));
         if (my_cut.colname() == "TIME") {
            if (my_cut.intervalType() == RangeCut::CLOSED) {
               gti = gti.applyTimeRangeCut(my_cut.minVal(), my_cut.maxVal());
            } else if (my_cut.intervalType() == RangeCut::MINONLY) {
               gti = gti.applyTimeRangeCut(my_cut.minVal(), gti.maxValue());
            } else if (my_cut.intervalType() == RangeCut::MAXONLY) {
               gti = gti.applyTimeRangeCut(gti.minValue(), my_cut.maxVal());
            }
         }
      }
   }
   gti.writeExtension(eventFile);
}

std::string CutController::filterString() const {
   std::string filter(m_cuts.filterString());
   if (filter != "") {
      filter += " && ";
   }
   filter += "gtifilter()";
   if (m_evclsFilter != "") {
      filter = m_evclsFilter + " && " + filter;
   }
   return filter;
}

} // namespace dataSubselector
