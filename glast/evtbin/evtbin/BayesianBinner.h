/** \file BayesianBinner.h
    \brief Declaration of BayesianBinner class.
*/
#ifndef evtbin_BayesianBinner_h
#define evtbin_BayesianBinner_h

#include <deque>
#include <string>
#include <vector>

#include "evtbin/OrderedBinner.h"

namespace evtbin {
  /** \class BayesianBinner
      \brief Declaration of a linearly uniform interval binner.
  */
  class BayesianBinner : public OrderedBinner {
    public:
      /** \brief Construct a Bayesian block binner object, using any kind of iterator to provide the cell population.
      */
      template <typename Itor>
      BayesianBinner(const IntervalCont_t & intervals, Itor cell_begin, const std::string & name = std::string(),
        double ncp_prior = 9.): OrderedBinner(IntervalCont_t(), name), m_cell_pop(cell_begin, cell_begin + intervals.size()),
        m_ncp_prior(ncp_prior) {
        computeBlocks(intervals);
      }

      virtual ~BayesianBinner() throw();

      /** \brief Create copy of this object.
      */
      virtual Binner * clone() const;

    private:
      /** \brief Perform the Bayesian Block procedure to determine the block definitions.
      */
      void computeBlocks(const IntervalCont_t & intervals);

      /** \brief Internal utility to compute the log posterior (Bayes factor).
      */
      void computeLogProb(const std::deque<double> & rev_csize, const std::deque<double> & rev_cpop,
        std::vector<double> & result) const;

      std::vector<double> m_cell_pop;
      const double m_ncp_prior;
  };

}

#endif
