/** \file Hist1D.h
    \brief One dimensional histogram.
*/
#ifndef evtbin_Hist1D_h
#define evtbin_Hist1D_h

#include <vector>

#include "evtbin/Hist.h"

namespace evtbin {

  class Binner;

  /** \class Hist1D
      \brief One dimensional histogram.
  */
  class Hist1D : public Hist {
    public:
      typedef std::vector<double> Cont_t;
      typedef Cont_t::const_iterator ConstIterator;

      /** \brief Create a one dimensional histogram which uses the given binner object:
          \param binner The binner object to use when filling bins.
      */
      Hist1D(const Binner & binner);

      virtual ~Hist1D() throw();

      /** \brief Increment the bin appropriate for the given value.
                 This is generic for N-dimensional histograms.
          \param value Vector giving the value being binned. The vector must have at least as
                 many values as the dimensionality of the histogram.
      */
      virtual void fillBin(const std::vector<double> & value, double weight = 1.);

      /** \brief Increment the bin appropriate for the given value.
          \param value The value being binned.
      */
      void fillBin(double value, double weight = 1.);

      const double & operator [](Cont_t::size_type index) const;

      ConstIterator begin() const;

      ConstIterator end() const;

    private:
      Cont_t m_data;
  };

  inline const double & Hist1D::operator [](Cont_t::size_type index) const { return m_data[index]; }

  inline Hist1D::ConstIterator Hist1D::begin() const { return m_data.begin(); }

  inline Hist1D::ConstIterator Hist1D::end() const { return m_data.end(); }

}

#endif
