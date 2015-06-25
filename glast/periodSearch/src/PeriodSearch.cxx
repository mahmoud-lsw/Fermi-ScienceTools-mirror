/** \file PeriodSearch.cxx
    \brief Implementation of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "PeriodSearch.h"

#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "StatisticViewer.h"

PeriodSearch::PeriodSearch(size_type num_bins, const std::string & freq_unit): m_freq_unit(freq_unit), m_description(),
  m_viewer(2, num_bins) {
  // Set default labels to the viewer.
  m_viewer.setLabel(0, "FREQUENCY");
  m_viewer.setLabel(1, "STATISTIC");

  // Set default unit to the viewer.
  m_viewer.setUnit(0, m_freq_unit);
}

void PeriodSearch::updateViewer(double min_freq, double max_freq) {
  // Find position of the maximum in the range.
  std::pair<double, double> max = findMax(min_freq, max_freq);

  // Compute probability for one trial.
  std::pair<double, double> chance_prob = computeChanceProbOneTrial(max.second);

  // Compute the number of independent trials.
  size_type num_indep_trials = computeNumIndepTrials(min_freq, max_freq);

  // Compute the multi-trial chance probability.
  chance_prob.first = computeChanceProbMultiTrial(chance_prob.first, num_indep_trials);
  chance_prob.second = computeChanceProbMultiTrial(chance_prob.second, num_indep_trials);

  // Compute number of bins.
  std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
  size_type num_bins = indices.second - indices.first;

  // Reset min/max frequency if either bound was not explicitly specified (negative).
  const StatisticViewer::data_type & freq = m_viewer.getData(0);
  if (0. > min_freq) min_freq = freq[indices.first];
  if (0. > max_freq && indices.second > 0) max_freq = freq[indices.second - 1];

  // Create caption for the viewer.
  std::ostringstream os;
  std::streamsize orig_precision = os.precision();
  os.precision(std::numeric_limits<double>::digits10);
  os << getDescription() << std::endl
     << "Search Range (" << m_freq_unit << "): [" << min_freq << ", " << max_freq << "]" << std::endl
     << "Number of Trial Frequencies: " << num_bins << std::endl
     << "Number of Independent Trials: " << num_indep_trials << std::endl
     << "Maximum Statistic: " << max.second << " at " << max.first << " " << m_freq_unit << std::endl
     << "Chance Probability Range: " << "(" << chance_prob.first << ", " << chance_prob.second << ")";
  os.precision(orig_precision);

  // Set caption to the viewer.
  m_viewer.setCaption(os.str());

  // Impose range limits.
  m_viewer.selectData(indices.first, indices.second);
}

std::pair<double, double> PeriodSearch::findMax(double min_freq, double max_freq) const {
  bool found_max = false;
  size_type max_idx = 0;
  double max = 0.;

  // Impose range limits.
  std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
  size_type begin_index = indices.first;
  size_type end_index = indices.second;

  for (size_type ii = begin_index; ii < end_index; ++ii) {
    // If the value is larger than the current maximum, replace it.
    const StatisticViewer::data_type & spec = m_viewer.getData(1);
    if (spec[ii] > max) {
      max = spec[ii];
      max_idx = ii;
      found_max = true;
    }
  }

  // Make sure a valid maximum was found.
  if (!found_max) {
    std::ostringstream os;
    os << "Could not find the maximum statistic at any trial frequency in range [" << min_freq << ", " << max_freq << "]";
    throw std::runtime_error(os.str());
  }
  const StatisticViewer::data_type & freq = m_viewer.getData(0);
  return std::pair<double, double>(freq[max_idx], max);
}

st_stream::OStream & PeriodSearch::write(st_stream::OStream & os) const {
  return os << getDescription();
}

double PeriodSearch::computeChanceProbMultiTrial(double prob_one_trial, size_type num_indep_trial) {
  static double epsilon = std::numeric_limits<double>::epsilon();
  double chance_prob = 0.;

  // Handle special cases.
  if (0 == num_indep_trial) chance_prob = 0.;
  else if (1 == num_indep_trial) chance_prob = prob_one_trial;
  else if (1. == prob_one_trial) chance_prob = 1.;
  else if (0. == prob_one_trial) chance_prob = 0.;
  else {
    // Compute x == -N * ln(1 - p) where N == num_indep_trial and p == prob_one_trial.
    double xx = 0.;
    if (.1 <= prob_one_trial) {
      // For large p use the formula directly.
      xx = -(signed long)(num_indep_trial) * std::log(1. - prob_one_trial);
    } else {
      // For small values of p, use iterative power series expansion.
      // Helper: p^n.
      double p_to_the_n = prob_one_trial;
      // qn == p^n/n -> q1 = p / 1.
      double q_sub_n = p_to_the_n;
      // Sn == sum qn -> Initial sum is 0.
      double S_sub_n = 0.;
      // Iterate over nn starting with 1.
      size_type nn = 1;
      do {
        // Compute Sn.
        S_sub_n += q_sub_n;
        // Prepare the next qn term and helper p_to_the_n:
        // n becomes "n+1"
        ++nn;
        // Compute p^"n+1"
        p_to_the_n *= prob_one_trial;
        // Compute q_sub_"n+1"
        q_sub_n = p_to_the_n / nn;
      } while (epsilon < q_sub_n / S_sub_n); // Note that this compares q_sub_n for the current nn to S_sub_n for the previous nn.

      // x == N * Sn for last n.
      xx = num_indep_trial * S_sub_n;
    }

    // Compute chance_prob == Q == 1. - exp(-x).
    if (.1 <= xx) {
      // For large x use the formula directly.
      chance_prob = 1. - std::exp(-xx);
    } else {
      // For small values of x, use iterative power series expansion.
      // Helper: leading term is x^(2n+1) / (2n + 1)!. For n = 0 this is x^1/1 = x.
      double leading_term = xx;
      // qn == leading term * (1 - xx/(2n + 2)). For n = 0 this is leading term * (1 - x/2)
      double q_sub_n = leading_term * (1. - xx / 2.);
      // Sn == sum qn -> Initial sum is 0.
      double S_sub_n = 0.;
      // Iterate over nn starting with 1.
      size_type nn = 0;
      do {
        // Compute Sn.
        S_sub_n += q_sub_n;
        // Prepare the next helper leading_term and qn.
        // n becomes "n+1"
        ++nn;
        // Compute next leading term.
        leading_term *= xx * xx / (2 * nn * (2 * nn + 1));
        // Compute q_sub_"n+1"
        q_sub_n = leading_term * (1. - xx / (2 * nn + 2));
      } while (epsilon < q_sub_n / S_sub_n); // Note that this compares q_sub_n for the current nn to S_sub_n for the previous nn.

      // Q == Sn for last n.
      chance_prob = S_sub_n;
    }
  }

  return chance_prob;
}

std::pair<PeriodSearch::size_type, PeriodSearch::size_type> PeriodSearch::getRangeIndex(double min_freq,
  double max_freq) const {

  const StatisticViewer::data_type & freq = m_viewer.getData(0);
  size_type begin_index = 0;
  size_type end_index = freq.size();

  if (0. <= min_freq) {
    // Find first element whose frequency is not less than the minimum frequency.
    for (begin_index = 0; begin_index < freq.size() && freq[begin_index] < min_freq; ++begin_index);
  }

  if (0. <= max_freq) {
    // Find last element whose frequency is not greater than the maximum frequency.
    for (end_index = freq.size(); end_index > begin_index && freq[end_index - 1] > max_freq; --end_index);
  }

  if (begin_index >= end_index) {
    std::ostringstream os;
    os << "No bins contained in frequency search range [" << min_freq << ", " << max_freq << "] " << m_freq_unit;
    throw std::runtime_error(os.str());
  }

  return std::make_pair(begin_index, end_index);
}

std::string PeriodSearch::getDescription() const {
  return m_description;
}

StatisticViewer & PeriodSearch::getViewer() {
  return m_viewer;
}

const StatisticViewer & PeriodSearch::getViewer() const {
  return m_viewer;
}

void PeriodSearch::setDescription(const std::string & search_type, double fourier_res, double sampling_step,
  const std::string & search_info) {
  // Write out common parameters.
  std::ostringstream os;
  os << "Search Type: " << search_type << std::endl
     << "Fourier Resolution: " << fourier_res << " " << m_freq_unit << std::endl
     << "Sampling Frequency: " << sampling_step << " " << m_freq_unit << std::endl
     << search_info;

  // Set the description.
  m_description = os.str();
}

st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test) {
  return test.write(os);
}
