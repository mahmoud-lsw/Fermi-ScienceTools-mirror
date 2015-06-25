/**
 * @file PowerLaw.h
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/src/test/PowerLaw.h,v 1.2 2004/08/26 17:02:24 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>

class PowerLaw {

public:

   PowerLaw(double prefactor, double index) : m_prefactor(prefactor),
                                              m_index(index) {}

   double operator()(double x) const {
      return m_prefactor*std::pow(x, m_index);
   }

   double index() const {
      return m_index;
   }

   double prefactor() const {
      return m_prefactor;
   }

private:

   double m_prefactor;
   double m_index;

};
