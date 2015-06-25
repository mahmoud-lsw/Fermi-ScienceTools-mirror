/** 
 * @file Bilinear.h
 * @brief Definition of interpolating class Bilinear.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/st_facilities/Bilinear.h,v 1.1 2012/12/10 10:17:45 cohen Exp $
 */

#ifndef st_facilities_Bilinear_h
#define st_facilities_Bilinear_h

#include <vector>

namespace st_facilities {

/**  
 * @class Bilinear
 *
 * @brief Bilinear interpolator.  
 *
 */

class Bilinear {

public:

   Bilinear(const std::vector<double> & x,
            const std::vector<double> & y, 
            const std::vector<double> & values, 
            double xlo, double xhi, double ylo, double yhi);

   double operator()(double x, double y) const;

//    void getCorners(double x, double y, 
//                    double & tt, double & uu,
//                    std::vector<double> & corner_xvals,
//                    std::vector<double> & corner_yvals,
//                    std::vector<double> & zvals) const;

//    static double evaluate(double tt, double uu, 
//                           const std::vector<double> & yvals);

   void getCorners(double x, double y, 
                   double & tt, double & uu,
                   double * corner_xvals,
                   double * corner_yvals,
                   double * zvals) const;

   static double evaluate(double tt, double uu, 
                          const double * zvals);

   double getPar(size_t i, size_t j) const;

   void setPar(size_t i, size_t j, double value);

private:

   std::vector<double> m_x;
   std::vector<double> m_y;
   std::vector<double> m_values;
   
};

} // namespace st_facilities

#endif // st_facilities_Bilinear_h
