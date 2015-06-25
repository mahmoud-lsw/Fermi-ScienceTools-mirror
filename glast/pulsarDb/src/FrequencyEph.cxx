/** \file FrequencyEph.cxx
    \brief Implementation of the FrequencyEph class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/FrequencyEph.h"

#include "timeSystem/CalendarFormat.h"
#include "timeSystem/MjdFormat.h"
#include "timeSystem/SourcePosition.h"

using namespace timeSystem;

namespace tip {
  class Header;
}

namespace pulsarDb {

  FrequencyEph::FrequencyEph(const tip::Table::ConstRecord & record, const tip::Header & /* header */):
    m_system(&TimeSystem::getSystem("TDB")), m_since("TDB", 0, 0.), m_until("TDB", 0, 0.), m_epoch("TDB", 0, 0.),
    m_ra(0.), m_dec(0.), m_phi0(0.), m_f0(0.), m_f1(0.), m_f2(0.) {
    // Epoch and toa are split into int and frac parts.
    long epoch_int = 0;
    double epoch_frac = 0.;
    long toa_int = 0;
    double toa_frac = 0.;

    // Read reference epoch and pulse TOA (integer parts required, fractional parts optional).
    read(record, "EPOCH_INT", epoch_int);
    read(record, "EPOCH_FRAC", epoch_frac, 0.);
    read(record, "TOABARY_INT", toa_int);
    read(record, "TOABARY_FRAC", toa_frac, 0.);

    // Combine separate parts of epoch and toa to get single values.
    m_epoch = AbsoluteTime("TDB", Mjd(epoch_int, epoch_frac));
    AbsoluteTime toa("TDB", Mjd(toa_int, toa_frac));

    // Read the start time of validity window (required).
    long valid_since_date = 0;
    read(record, "VALID_SINCE", valid_since_date);
    m_since = AbsoluteTime("TDB", Mjd(valid_since_date, 0.));

    // Read the end time of validity window (required).
    // Note: One is added to the endpoint because the "VALID_UNTIL" field in the file expires at the end of that day,
    // whereas the valid_until argument to the ephemeris object is the absolute cutoff.
    long valid_until_date = 0;
    read(record, "VALID_UNTIL", valid_until_date);
    m_until = AbsoluteTime("TDB", Mjd(valid_until_date + 1, 0.));

    // Read the sky position and frequency coefficients (RA, DEC, F0: required, F1, F2: optional).
    read(record, "RA",  m_ra);
    read(record, "DEC", m_dec);
    read(record, "F0",  m_f0);
    read(record, "F1",  m_f1 , 0.);
    read(record, "F2",  m_f2 , 0.);

    // Create temporary copy of this ephemeris with phi0 == 0.
    FrequencyEph tmp("TDB", m_since, m_until, m_epoch, m_ra, m_dec, 0., m_f0, m_f1, m_f2);

    // Use the temporary ephemeris to compute the phase from the negative of the toa field.
    m_phi0 = - tmp.calcPulsePhase(toa);

    // Make sure it is in the range [0, 1). calcPulsePhase is bounded in this way.
    if (0. > m_phi0) m_phi0 += 1.;
  }

  void FrequencyEph::writeModelParameter(st_stream::OStream & os) const {
    std::string epoch_string;
    try {
      epoch_string = m_epoch.represent(m_system->getName(), MjdFmt);
    } catch (const std::exception &) {
      epoch_string = m_epoch.represent(m_system->getName(), CalendarFmt);
    } 
    os << format("Epoch", epoch_string, "")        << std::endl;
    os << format("RA",    m_ra,         "degrees") << std::endl;
    os << format("Dec",   m_dec,        "degrees") << std::endl;
    os << format("Phi0",  m_phi0,       "")        << std::endl;
    os << format("F0",    m_f0,         "s**(-1)") << std::endl;
    os << format("F1",    m_f1,         "s**(-2)") << std::endl;
    os << format("F2",    m_f2,         "s**(-3)");
  }

  double FrequencyEph::calcPulsePhase(const AbsoluteTime & ev_time, double phase_offset) const {
    // Compute a pulse phase value.
    double dt = calcElapsedSecond(ev_time);
    double dt_squared = dt * dt;
    double phase = m_phi0 + m_f0 * dt + m_f1/2.0 * dt_squared + m_f2/6.0 * dt * dt_squared;

    // Express phase as a value between 0. and 1., after adding a global phase offset.
    return trimPhaseValue(phase, phase_offset);
  }

  double FrequencyEph::calcFrequency(const AbsoluteTime & ev_time, int derivative_order) const {
    double return_value = 0.;
    if (0 > derivative_order) {
      std::ostringstream os;
      os << "Negative order of frequency derivative is given: " << derivative_order;
      throw std::runtime_error(os.str());
    } else if (0 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = m_f0 + m_f1 * dt + 0.5 * m_f2 * dt * dt;
    } else if (1 == derivative_order) {
      double dt = calcElapsedSecond(ev_time);
      return_value = m_f1 + m_f2 * dt;
    } else if (2 == derivative_order) {
      return_value = m_f2;
    } else {
      return_value = 0.;
    }
    return return_value;
  }

  SourcePosition FrequencyEph::calcPosition(const AbsoluteTime & /* ev_time */) const {
    return SourcePosition(m_ra, m_dec);
  }

}
