/** 
 * @file BrokenPowerLaw.cxx
 * @brief Implementation for the BrokenPowerLaw Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/BrokenPowerLaw.cxx,v 1.6 2015/03/03 18:03:49 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "BrokenPowerLaw.h"

namespace optimizers {

BrokenPowerLaw::
BrokenPowerLaw(double Prefactor, double Index1,
               double Index2, double BreakValue) 
   : Function("BrokenPowerLaw", 4, "Prefactor", "dArg", Addend) {
   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index1"), Index1, true);
   addParam(std::string("Index2"), Index2, true);
   addParam(std::string("BreakValue"), BreakValue, true);
}

double BrokenPowerLaw::value(Arg & xarg) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index1, Index2, BreakValue};

   std::vector<Parameter> my_params;
   getParams(my_params);

   if (x < my_params[BreakValue].getTrueValue()) {
      return my_params[Prefactor].getTrueValue()
         *pow((x/my_params[BreakValue].getTrueValue()), 
              my_params[Index1].getTrueValue());
   } else {
      return my_params[Prefactor].getTrueValue()
         *pow((x/my_params[BreakValue].getTrueValue()), 
              my_params[Index2].getTrueValue());
   }
   return 0;
}

double BrokenPowerLaw::derivByParamImp(Arg & xarg, 
                                       const std::string & paramName) const {

   double x = dynamic_cast<dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index1, Index2, BreakValue};

   std::vector<Parameter> my_params;
   getParams(my_params);

   double gam1 = -my_params[Index1].getTrueValue();
   double gam2 = -my_params[Index2].getTrueValue();
   double Eb = my_params[BreakValue].getTrueValue();

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
      throw ParameterNotFound(paramName, 
                              getName(), 
                              "BrokenPowerLaw::derivByParam");
   }
   
   switch (iparam) {
   case Prefactor:
      if (x < Eb) {
         return std::pow(x/Eb, -gam1)
            *my_params[Prefactor].getScale();
      } else {
         return std::pow(x/Eb, -gam2)
            *my_params[Prefactor].getScale();
      }
      break;
   case Index1:
      if (x < Eb) {
         return value(xarg)*std::log(x/Eb)
            *my_params[Index1].getScale();
      } else {
         return 0;
      }
      break;
   case Index2:
      if (x < my_params[BreakValue].getTrueValue()) {
         return 0;
      } else {
         return value(xarg)*std::log(x/Eb)
            *my_params[Index2].getScale();
      }
      break;
   case BreakValue:
      if (x < Eb) {
         return value(xarg)*gam1/Eb
            *my_params[BreakValue].getScale();
      } else {
         return value(xarg)*gam2/Eb
            *my_params[BreakValue].getScale();
      }
      break;
   default:
      break;
   }
   return 0;
}

} // namespace optimizers
