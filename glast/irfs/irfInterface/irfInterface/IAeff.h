/** 
 * @file IAeff.h
 * @brief IAeff class definition.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/irfInterface/IAeff.h,v 1.6 2009/05/30 22:28:22 jchiang Exp $
 */

#ifndef irfInterface_IAeff_h
#define irfInterface_IAeff_h

namespace astro {
   class SkyDir;
}

namespace irfInterface {

/** 
 * @class IAeff
 *
 * @brief Abstract interface for the LAT effective area classes.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfInterface/irfInterface/IAeff.h,v 1.6 2009/05/30 22:28:22 jchiang Exp $
 */

class IAeff {
    
public:

   virtual ~IAeff() {}

   /// Pure virtual method to define the interface for the member
   /// function returning the effective area in cm^2.
   /// @param energy True photon energy (MeV).
   /// @param srcDir True photon direction.
   /// @param scZAxis Spacecraft z-axis.
   /// @param scXAxis Spacecraft x-axis.
   /// @param time Photon arrival time (MET s).
   virtual double value(double energy, 
                        const astro::SkyDir &srcDir, 
                        const astro::SkyDir &scZAxis,
                        const astro::SkyDir &scXAxis,
                        double time=0) const = 0;

   /// Return the effective area (cm^2) as a function of instrument 
   /// coordinates.
   /// @param energy True photon energy (MeV).
   /// @param theta True inclination angle (degrees).
   /// @param phi True azimuthal angle measured wrt the instrument
   ///             X-axis (degrees).
   /// @param time Photon arrival time (MET s).
   virtual double value(double energy, double theta, double phi,
                        double time=0) const = 0;

   /// This method is also virtual, in case the sub-classes wish to
   /// overload it.
   virtual double operator()(double energy, 
                             const astro::SkyDir &srcDir, 
                             const astro::SkyDir &scZAxis,
                             const astro::SkyDir &scXAxis,
                             double time=0) const {
      return value(energy, srcDir, scZAxis, scXAxis, time);}

   virtual IAeff * clone() = 0;

   /// @return An absolute upper limit on the value of the effective
   /// area for all energies and directions (cm^2).
   virtual double upperLimit() const = 0;

   virtual void setPhiDependence(bool usePhiDependence) {
      m_usePhiDependence = usePhiDependence;
   }

   virtual bool usePhiDependence() const {
      return m_usePhiDependence;
   }

protected:

   bool m_usePhiDependence;

};

} // namespace irfInterface

#endif // irfInterface_IAeff_h
