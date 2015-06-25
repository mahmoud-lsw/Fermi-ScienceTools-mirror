/**
 * @file WcsMap2.h
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/WcsMap2.h,v 1.1.1.5.2.5 2015/04/26 06:53:44 jasercio Exp $
 */

#ifndef Likelihood_WcsMap2_h
#define Likelihood_WcsMap2_h

#include <vector>
#include <utility>

#include "astro/SkyDir.h"

namespace Likelihood {

class BinnedExposure;
class DiffuseSource;
class MeanPsf;

/**
 * @class WcsMap2
 *
 */

class WcsMap2 {

public:

   WcsMap2(const std::string & filename, const std::string & extension="",
           bool interpolate=true, bool enforceEnergyRange=false);

   WcsMap2(const DiffuseSource & diffuseSource, double ra, double dec,
           double pix_size, int npts, double energy=100.,
           const std::string & proj_name="STG", bool use_lb=false,
           bool interpolate=false, bool enforceEnergyRange=false);

   WcsMap2(const DiffuseSource & diffuseSource, double ra, double dec,
           double crpix1, double crpix2, double cdelt1, double cdelt2,
           int naxis1, int naxis2, double energy=100.,
           const std::string & proj_name="STG", bool use_lb=false,
           bool interpolate=false, bool enforceEnergyRange=false);

   ~WcsMap2();

   WcsMap2(const WcsMap2 &);

   WcsMap2 & operator=(const WcsMap2 &);

   double operator()(const astro::SkyDir & dir, double energy=-1) const;

   double operator()(const astro::SkyDir & dir, int k) const;

   WcsMap2 convolve(double energy, const MeanPsf & psf,
                    const BinnedExposure & exposure,
                    bool performConvolution=true,
                    int k=0) const;

//    const std::vector< std::vector< std::vector<double> > > & image() const {
//       return m_image;
//    }
   const std::vector< std::vector< std::vector<float> > > & image() const {
      return m_image;
   }


   /// @return Solid angle of the (ilon, ilat) pixel
   static double solidAngle(const astro::SkyProj & proj, 
                            double ilon, double ilat);

   double solidAngle(double ilon, double ilat) const;

//   const std::vector< std::vector<double> > & solidAngles() const;
   const std::vector< std::vector<float> > & solidAngles() const;

   /// @return Pixel value as a function of index
   double pixelValue(double ilon, double ilat, int k=0) const;
   
   /// @return SkyDir corresponding to the pixel indices
   astro::SkyDir skyDir(double ilon, double ilat) const;

   int nxpix() const {
      return m_naxis1;
   }

   int nypix() const {
      return m_naxis2;
   }

   int nenergies() const {
      return m_naxis3;
   }

   const std::vector<double> & energies() const {
      return m_energies;
   }

   bool insideMap(const astro::SkyDir & dir) const;

   bool withinMapRadius(const astro::SkyDir & dir) const;

   std::pair<astro::SkyDir, astro::SkyDir> 
   minMaxDistPixels(const astro::SkyDir & dir) const;

   void getCorners(std::vector<astro::SkyDir> & corners) const;

   double mapIntegral() const;

   double mapIntegral(double energy) const;

   double cdelt1() const {
      return m_cdelt1;
   }

   double cdelt2() const {
      return m_cdelt2;
   }

   /// Rebin the map data by the specified factor.  By default, this
   /// will average over combined pixels, as would be appropriate for
   /// intensity or exposure maps.  For counts maps, set
   /// average=false, so that combined bins are summed.  Return a new
   /// WcsMap2 object with the new geometry.  The reference direction
   /// will be unchanged from the original (and so will not generally
   /// point to the center of the map.)
   WcsMap2 * rebin(unsigned int factor, bool average=true);

   void setExtrapolation(bool enforceEnergyRange);

   bool enforceEnergyRange() const;

   unsigned long extrapolated() const {
      return m_extrapolated;
   }

   static double interpolatePowerLaw(double x, double x1, double x2,
                                     double y1, double y2);

private:

   astro::SkyDir m_refDir;

//   typedef std::vector< std::vector<double> > ImagePlane_t;
   typedef std::vector< std::vector<float> > ImagePlane_t;

   std::vector<ImagePlane_t> m_image;

   mutable ImagePlane_t m_solidAngles;

   int m_naxes;
   int m_naxis1;
   int m_naxis2;
   int m_naxis3;

   astro::SkyProj * m_proj;

   /// astro::SkyProj provides almost no introspection, so we store the
   /// map projection locally.
   double m_crpix1, m_crpix2;
   double m_crval1, m_crval2;
   double m_cdelt1, m_cdelt2;
   double m_crota2;
   
   std::vector<double> m_energies;

   bool m_interpolate;

   bool m_isPeriodic;

   astro::SkyDir::CoordSystem m_coordSys;

   double m_mapIntegral;

//   std::vector<double> m_mapIntegrals;
   std::vector<float> m_mapIntegrals;

   bool m_enforceEnergyRange;

   mutable unsigned long m_extrapolated;

   std::string m_filename;

   // enclosing cone half angle for map
   double m_mapRadius;

   WcsMap2();

   void computeMapIntegrals();

   void check_energy_index(int k) const;

   void check_energy(double energy) const;

   void check_negative_pixels(const ImagePlane_t &) const;

};

} // namespace Likelihood

#endif // Likelihood_WcsMap2_h