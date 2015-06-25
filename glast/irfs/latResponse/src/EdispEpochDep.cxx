/**
 * @file EdispEpochDep.cxx
 * @brief Definition for class to handle epoch-dependent energy dispersion.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * 
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/EdispEpochDep.cxx,v 1.5.2.4 2015/04/26 06:53:49 jasercio Exp $
 */

#include "EdispEpochDep.h"

namespace latResponse {

EdispEpochDep::EdispEpochDep() {
}

EdispEpochDep::EdispEpochDep(const EdispEpochDep & other) 
   : irfInterface::IEdisp(other), EpochDep(other) {
   for (size_t i(0); i < other.m_edisps.size(); i++) {
      m_edisps.push_back(other.m_edisps[i]->clone());
   }
}

EdispEpochDep & EdispEpochDep::operator=(const EdispEpochDep & rhs) {
   irfInterface::IEdisp::operator=(rhs);
   EpochDep::operator=(rhs);
   for (size_t i(0); i < m_edisps.size(); i++) {
      delete m_edisps[i];
   }
   m_edisps.clear();
   for (size_t i(0); i < rhs.m_edisps.size(); i++) {
      m_edisps.push_back(rhs.m_edisps[i]->clone());
   }
}

EdispEpochDep::~EdispEpochDep() {
   for (size_t i(0); i < m_edisps.size(); i++) {
      delete m_edisps[i];
   }
}

double EdispEpochDep::value(double appEnergy, 
                            double energy,
                            const astro::SkyDir & srcDir, 
                            const astro::SkyDir & scZAxis,
                            const astro::SkyDir & scXAxis,
                            double time) const {
   size_t indx(index(time));
   return m_edisps[indx]->value(appEnergy, energy, srcDir, scZAxis,
                                scXAxis, time);
}

double EdispEpochDep::value(double appEnergy, double energy,
                            double theta, double phi,
                            double time) const {
   size_t indx(index(time));
   return m_edisps[indx]->value(appEnergy, energy, theta, phi, time);
}

void EdispEpochDep::addEdisp(const irfInterface::IEdisp & edisp,
                             double epoch_start) {
   appendEpoch(epoch_start);
   m_edisps.push_back(const_cast<irfInterface::IEdisp &>(edisp).clone());
}

} // namespace latResponse
