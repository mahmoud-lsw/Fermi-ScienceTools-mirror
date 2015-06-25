/** \file ChiSquaredProb.cxx
    \brief Implmementation for ChiSquaredProb class.
    \author Masaharu Hirayama, GSSC
*/
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "ChiSquaredProb.h"

ChiSquaredProb::ChiSquaredProb(int dof, double min_pdf): m_dof(dof), m_dof_minus_2(dof-2), m_lognorm(0.0) {
  // Check the sign of dof argument.
  if (dof <= 0) {
    std::ostringstream os;
    os << "Non-positive number is given for the degrees of freedom of a chi-squared distribution: " << dof;
    throw std::logic_error(os.str());
  }

  // Check the sign of min_pdf argument.
  if (min_pdf <= 0.) {
    std::ostringstream os;
    os << "Non-positive number is given for the minimum value of probability density function to compute: " << min_pdf;
    throw std::logic_error(os.str());
  }

  // Pre-compute normalization.
  double nterm;
  if (dof % 2) {
    m_lognorm = - 0.5 * std::log(M_PI);
    nterm = dof / 2;
  } else {
    m_lognorm = 0.0;
    nterm = dof / 2 - 2;
  }
  for (int ii=1; ii<=nterm; ii++) {
    m_lognorm -= std::log(m_dof / 2.0 - static_cast<double>(ii));
  }
  m_lognorm -= m_dof / 2.0 * M_LN2;

  // Find the maximum chi-squared value to compute.
  // NOTE: in the following search for m_max_chisq, the initial value
  // of x_cur must be > m_dof_minus_2 and > 0.
  double x_cur = m_dof + 1.0;
  if (pdf(x_cur) > min_pdf) {
    double log_min_pdf = (std::log(min_pdf) - m_lognorm) * 2.0;
    const double MY_EPSILON = 1.0e-10;
    const double min_diff = 2.0 * std::log(1.0 + MY_EPSILON);
    double log_diff;
    while (1) {
      log_diff = m_dof_minus_2 * std::log(x_cur) - x_cur - log_min_pdf;
      if (fabs(log_diff) > min_diff) x_cur += log_diff;
      else break;
    }
  }
  m_max_chisq = x_cur;
}

std::pair<double,double> ChiSquaredProb::operator() (double chisq, double precision, int max_iteration) const {
  // Return 1.0 in trivial cases.
  if (chisq <= 0.0) return std::make_pair(1.0, 1.0);

  // Return the maximum residual if chisq >= m_max_chisq.
  if (chisq >= m_max_chisq) {
    if (m_max_chisq > m_dof) {
      double max_residual = 2.0 * m_max_chisq * pdf(m_max_chisq)
	/ (m_max_chisq - m_dof);
      return std::make_pair(0.0, (1.0 < max_residual ? 1.0 : max_residual));
    } else {
      return std::make_pair(0.0, 1.0);
    }
  }

  // Prepare for computing step size.
  double p_value = precision / 2.0;
  double x_offset = 0.0;
  if (m_dof > 2.1) {
    x_offset = 2.0 * p_value + std::sqrt(2.0 * p_value * m_dof_minus_2);
  }

  // Set initial values for numerical integration.
  double x_cur = chisq;
  double f_cur = pdf(x_cur);
  double Qmin = 0.0;
  double Qmax = 0.0;
  double Rn = 1.0;

  // Prepare Helper variables for numerical integration.
  double x_step, x_next, f_next;
  int num_iter = 0;

  // Integration where pdf(x) increases with x.
  if (m_dof > 2.1) {
    for (; (num_iter < max_iteration) && (x_cur < m_dof_minus_2)
	   && (x_cur < m_max_chisq);
	 num_iter++, x_cur = x_next, f_cur = f_next) {

      // Set next sampling point.
      x_step = p_value * fabs(eta(x_cur));
      x_next = x_cur + x_step;
      if (x_next > m_dof_minus_2) {
	x_next = m_dof_minus_2;
	x_step = x_next - x_cur;
      }
      f_next = pdf(x_next);

      // Accumulate estimators.
      Qmax += f_next * x_step;
      Qmin += f_cur  * x_step;
    }
  }

  // Integration where pdf(x) decreases with x.
  for (; (num_iter < max_iteration) && (x_cur < m_max_chisq);
       num_iter++, x_cur = x_next, f_cur = f_next) {
    // Set next sampling point.
    x_step = p_value * fabs(eta(x_cur + x_offset));
    x_next = x_cur + x_step;
    if (x_next > m_max_chisq) {
      x_next = m_max_chisq;
      x_step = x_next - x_cur;
    }
    f_next = pdf(x_next);

    // Accumulate estimators.
    Qmax += f_cur  * x_step;
    Qmin += f_next * x_step;

    if (x_next > m_dof) {
      // Compute residual.
      Rn = 2.0 * x_next * f_next / (x_next - m_dof);

      // Terminate numerical integration.
      if (Rn <= Qmin* p_value) break;
    }
  }

  // Add residual (upper limit).
  Qmax += Rn;
  Qmax = (Qmax < 1.0 ? Qmax : 1.0);

  // Return integration and upper limit for residual.
  return std::make_pair(Qmin, Qmax);
}
