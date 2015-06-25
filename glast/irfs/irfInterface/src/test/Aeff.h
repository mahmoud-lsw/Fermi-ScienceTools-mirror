/** 
 * @file Aeff.h
 * @brief Concrete IAeff subclass.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/src/test/Aeff.h,v 1.3 2006/08/14 19:34:40 jchiang Exp $
 */

#ifndef irfInterface_Aeff_h
#define irfInterface_Aeff_h

#include "irfInterface/IAeff.h"

/** 
 * @class Aeff
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/src/test/Aeff.h,v 1.3 2006/08/14 19:34:40 jchiang Exp $
 */

class Aeff : public irfInterface::IAeff {
    
public:

   Aeff() {}

   virtual ~Aeff() {}

   virtual double value(double, 
                        const astro::SkyDir &,
                        const astro::SkyDir &,
                        const astro::SkyDir &,
                        double) const {return 0;}

   virtual double value(double, double, double, double) const {return 0;}

   virtual Aeff * clone() {return new Aeff(*this);}

   virtual double upperLimit() const {return 0;}

};

#endif // irfInterface_Aeff_h
