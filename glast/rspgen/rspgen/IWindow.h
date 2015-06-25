/** \file IWindow.h
    \brief Abstract interface for window (region) specification.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_IWindow_h
#define rspgen_IWindow_h

#include <string>

#include "astro/SkyDir.h"

#include "irfInterface/IPsf.h"

namespace rspgen {

  /** \class IWindow
      \brief Abstract interface for window (region) specification.
  */
  class IWindow {
    public:
      virtual ~IWindow() throw() {}

      /// \brief Make a copy of this object.
      virtual IWindow * clone() const = 0;

      /** \brief Integrate the psf over this window for the given values.
          \param psf The psf object obtained from the irfs package.
          \param true_energy The true energy for which to evaluate the integral.
          \param theta The true inclination angle of the incident photon.
          \param phi The true azimuthal angle of the incident photon.
      */
      virtual double integrate(irfInterface::IPsf * psf, double true_energy, double theta, double phi) const = 0;
  };

}

#endif
