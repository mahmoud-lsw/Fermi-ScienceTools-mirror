/** 
 * @file FitsImage.h
 * @brief Declaration of FitsImage class
 * @authors J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/st_facilities/FitsImage.h,v 1.8 2006/01/18 05:47:38 jchiang Exp $
 *
 */

#ifndef st_facilities_FitsImage_h
#define st_facilities_FitsImage_h

#include <string>
#include <vector>

namespace astro {
   class SkyDir;
   class SkyProj;
}

namespace tip {
   class Header;
}

namespace st_facilities {

/** 
 * @class FitsImage
 *
 * @brief A class for accessing FITS image data.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/st_facilities/FitsImage.h,v 1.8 2006/01/18 05:47:38 jchiang Exp $
 *
 */

class FitsImage {
    
public:

   FitsImage() {}

   FitsImage(const std::string & fitsfile, 
             const std::string & extension="");

   FitsImage(const FitsImage &);

   virtual ~FitsImage();

   FitsImage & operator=(const FitsImage &);

   /// A vector of the image axes dimensions
   virtual void getAxisDims(std::vector<int> & axisDims) const;

   /// The names (CTYPEs) of the image axes
   virtual void getAxisNames(std::vector<std::string> & axisNames) const;

   /// This method computes arrays of longitude and latitude obtained
   /// by traversing the image plane by column number then row.
   virtual void getCelestialArrays(std::vector<double> & lonArray,
                                   std::vector<double> & latArray) const;
   
   /// Get the pixel values.  They will be indexed by column, row,
   /// then plane, i.e., indx = i + j*NAXIS1 + k*NAXIS1*NAXIS2.  Note
   /// that each image plane is indexed starting at the lower left
   /// (South-East) corner.
   virtual void getImageData(std::vector<double> & imageData) const;

   /// This returns the pixel solid angles.  Use of this method assumes
   /// that m_axis[0] represents a longitudinal coordinate and that
   /// m_axis[1] represents a latitudinal coordinate.  The pixel values
   /// will be indexed by column then row, indx = i + j*NAXIS1.
   virtual void getSolidAngles(std::vector<double> & solidAngles) const;

   const std::vector<float> & imageData() const {
      return m_image;
   }

   /// @brief Factory method to create an astro::SkyProj object.
   /// @param fitsFile FITS file containing the WCS projection information.
   /// @param extension The name of the relevent FITS image extension.
   static astro::SkyProj * skyProjCreate(const std::string & fitsFile,
                                         const std::string & extension="");

   /// @brief Compute the solid angle of a region defined by four
   /// ordered directions on the sky, joined by great circles.
   static double solidAngle(const astro::SkyDir & A,
                            const astro::SkyDir & B, 
                            const astro::SkyDir & C,
                            const astro::SkyDir & D);

protected:

/** 
 * @class AxisParams
 * @brief Nested n-tuple class to represent FITS image axis information
 */
   class AxisParams {
   public:
      AxisParams() {}
      ~AxisParams() {}
      int size;
      float refVal;
      float step;
      float refPixel;
      std::string axisType;
      std::string comment;
      bool logScale;
   };

   void read_fits_image();

   /// FITS file name.
   std::string m_filename;

   /// FITS extension name
   std::string m_extension;

   /// Descriptions for each image axis.
   std::vector<AxisParams> m_axes;

   /// The FITS image data
   std::vector<float> m_image;

   void setProjection(const tip::Header & header);

};

} // namespace st_facilities

#endif // st_facilities_FitsImage.h
