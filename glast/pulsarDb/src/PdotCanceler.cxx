/** \file PdotCanceler.cxx
    \brief Implementation of the PdotCanceler class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include "pulsarDb/PdotCanceler.h"

#include "pulsarDb/PulsarEph.h"

#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/TimeInterval.h"
#include "timeSystem/TimeSystem.h"

namespace pulsarDb {

  PdotCanceler::PdotCanceler(const std::string & time_system_name, const timeSystem::AbsoluteTime & time_origin,
    const std::vector<double> & fdot_ratio): m_time_system(&timeSystem::TimeSystem::getSystem(time_system_name)),
    m_time_origin(time_origin), m_fdot_ratio(fdot_ratio) {}

  PdotCanceler::PdotCanceler(const timeSystem::AbsoluteTime & time_origin, const PulsarEph & pulsar_eph, int max_derivative):
    m_time_system(&timeSystem::TimeSystem::getSystem(pulsar_eph.getSystem().getName())), m_time_origin(time_origin),
    m_fdot_ratio(max_derivative, 0.) {
    // Compute frequency derivatives.
    double f0 = pulsar_eph.calcFrequency(m_time_origin, 0);
    for (std::vector<double>::size_type ii = 0; ii < m_fdot_ratio.size(); ++ii) {
      m_fdot_ratio[ii] = pulsar_eph.calcFrequency(m_time_origin, ii + 1) / f0;
    }
  }

  const timeSystem::TimeSystem & PdotCanceler::getSystem() const {
    return *m_time_system;
  }

  void PdotCanceler::cancelPdot(timeSystem::AbsoluteTime & abs_time) const {
    // Compute elapsed seconds from time origin.
    const std::string time_system_name = m_time_system->getName();
    double dt = (abs_time - m_time_origin).computeDuration(time_system_name, "Sec");

    // Compute time correction.
    double correction = 0.;
    double factor = dt;
    for (std::vector<double>::size_type ii = 0; ii < m_fdot_ratio.size(); ++ii) {
      factor *= dt / (ii + 2);
      correction += factor * m_fdot_ratio[ii];
    }

    // Apply time correction.
    abs_time += timeSystem::ElapsedTime(time_system_name, timeSystem::Duration(correction, "Sec"));
  }

}
