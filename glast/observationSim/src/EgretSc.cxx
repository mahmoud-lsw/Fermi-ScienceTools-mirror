/**
 * @file EgretSc.cxx
 * @brief Implementation for EGRET spacecraft class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/observationSim/src/EgretSc.cxx,v 1.3 2010/06/03 04:16:37 jchiang Exp $
 */

#include <iostream>
#include <vector>
#include "EgretSc.h"

namespace observationSim {

CLHEP::HepRotation EgretSc::InstrumentToCelestial(double) {

// This implementation *should* ensure that an orthogonal set of axes
// are fed to the HepRotation constructor.
   CLHEP::Hep3Vector z_axis = m_zAxis();
   CLHEP::Hep3Vector x_axis = m_xAxis();
   CLHEP::Hep3Vector yAxis = z_axis.cross(x_axis);

   return CLHEP::HepRotation(yAxis.cross(z_axis), yAxis, z_axis);
}

} // namespace observationSim

