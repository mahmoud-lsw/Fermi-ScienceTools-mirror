// $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/astro/astro/PointingTransform.h,v 1.4 2006/03/21 01:43:17 usher Exp $
#ifndef OrbitModel_PointingTransform_H
#define OrbitModel_PointingTransform_H


// Include files
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "astro/SkyDir.h"
#include <cmath>
#include <utility> // for pair

namespace astro {


    /** @class PointingTransform
    * @brief Describe the GLAST pointing, and allow for transformations
    * @author S. Robinson 
    * Hold the GLAST pointing information, and allow for doing the transformation
    * between coordinate systems.  Generally, it is assumed that the user will want
    * to take a GLAST local incoming direction (x,y,z) and get the SkyDir direction
    * corresponding to that same direction.
    * <br>$Id: PointingTransform.h,v 1.4 2006/03/21 01:43:17 usher Exp $
    *
    */
    class PointingTransform
    {
    public:
        ///Constructors
        ///Z,X axis instantiation
        PointingTransform(SkyDir zdir, SkyDir xdir);

        /// The rotation that turns glast-local to SkyDir cartesian celestial vector
        CLHEP::HepRotation localToCelestial () const;

        /// The absolute direction corresponding to some GLAST direction
        SkyDir gDir(CLHEP::Hep3Vector localDir) const;

    private:
        SkyDir m_xDir,m_zDir;
    };


} // namespace astro
#endif    // LHCBEVENT_PointingTransform_H
