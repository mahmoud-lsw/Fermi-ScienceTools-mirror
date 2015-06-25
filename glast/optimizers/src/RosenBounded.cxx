/** 
 * @file RosenBounded.cxx
 * @brief Implementation for the 2D RosenBoundedbrock objective function
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/optimizers/src/RosenBounded.cxx,v 1.1.2.4 2015/04/26 06:53:54 jasercio Exp $
 */

#include <cmath>

#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"
#include "RosenBounded.h"

namespace optimizers {

RosenBounded::RosenBounded(double prefactor) 
   : Statistic("RosenBounded", 2), m_prefactor(prefactor), m_xmin(-1e30),
     m_xmax(1e30), m_ymin(-1e30), m_ymax(1e30) {
   addParam(std::string("x"), 1, true);
   addParam(std::string("y"), 1, true);
}

double RosenBounded::value(Arg &) const {
   double x = m_parameter[0].getTrueValue();
   double y = m_parameter[1].getTrueValue();

   check_bounds(x, y);

   return -(m_prefactor*pow((y - x*x), 2) + pow((1 - x), 2));
}

double RosenBounded::derivByParamImp(Arg &, 
                                     const std::string & paramName) const {
   std::vector<double> params;
   getParamValues(params);

   double x = params[0];
   double y = params[1];

   check_bounds(x, y);

   if (paramName == "x") {
      return -(-4.*m_prefactor*(y - x*x)*x - 2.*(1. - x));
   } else if (paramName == "y") {
      return -2.*m_prefactor*(y - x*x);
   }
   throw ParameterNotFound(paramName, getName(), "RosenBounded::derivByParam");
}

void RosenBounded::set_xbounds(double xmin, double xmax) {
   m_xmin = xmin;
   m_xmax = xmax;
}

void RosenBounded::set_ybounds(double ymin, double ymax) {
   m_ymin = ymin;
   m_ymax = ymax;
}

void RosenBounded::check_bounds(double x, double y) const {
   if (x < m_xmin || x > m_xmax) {
      std::ostringstream message;
      message << "RosenBound::value: "
              << "attempt to evaluate outside x-bounds. "
              << "x = " << x << std::endl;
      throw std::runtime_error(message.str());
   }
   if (y < m_ymin || y > m_ymax) {
      std::ostringstream message;
      message << "RosenBound::value: "
              << "attempt to evaluate outside y-bounds. "
              << "y = " << y << std::endl;
      throw std::runtime_error(message.str());
   }
}

} // namespace optimizers
