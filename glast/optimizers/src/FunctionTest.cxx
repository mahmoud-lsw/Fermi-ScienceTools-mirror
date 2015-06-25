/**
 * @file FunctionTest.cxx
 * @brief Unit test implementation for Function class hierarchy
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/FunctionTest.cxx,v 1.12 2013/01/10 16:31:48 jchiang Exp $
 */

#include <cassert>
#include <iostream>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "optimizers/FunctionTest.h"
#include "optimizers/Exception.h"

namespace {
   void do_nothing(double x) {
      x;
   }

   double select_delta(double x, double h) {
      double temp = x + h;
      do_nothing(temp);
      h = temp - x;
      return h;
   }
} // anonymous namespace

namespace optimizers {

void FunctionTest::parameters(const std::vector<Parameter> &params) {

   std::vector<std::string> paramNames;
   m_func->getParamNames(paramNames);
   unsigned int numParams = m_func->getNumParams();
   assert(numParams == paramNames.size());

// set the Parameter values
   for (unsigned int i = 0; i < numParams && i < params.size(); i++) {
      m_func->setParam(paramNames[i], params[i].getValue());
      m_func->parameter(paramNames[i]).setBounds(params[i].getBounds().first,
                                                 params[i].getBounds().second);
   }

// and test in groups and individually
   std::vector<double> paramValues;
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      assert(paramValues[i] == params[i].getValue());
      assert(m_func->getParamValue(paramNames[i]) == params[i].getValue());
   }

// set Parameters using Parameter class
   std::vector<Parameter> my_params;
   for (unsigned int i = 0; i < numParams && i < params.size(); i++) {
      my_params.push_back(params[i]);
// ensure names match up 
      my_params[i].setName(paramNames[i]);
   }

// set Parameters as a group
   m_func->setParams(my_params);

// test values again
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      assert(paramValues[i] == params[i].getValue());
      assert(m_func->getParamValue(paramNames[i]) == params[i].getValue());
   }

// reset Parameters individually with different values
   for (unsigned int i = 0; i < my_params.size(); i++) {
      my_params[i].setValue(params[i].getValue()/10.);
      m_func->setParam(my_params[i]);
   }

// test yet again
   m_func->getParamValues(paramValues);
   for (unsigned int i = 0; i < paramValues.size(); i++) {
      double value = params[i].getValue()/10.;
      if (paramValues[i] != value) {
         // std::cout << paramValues[i] << "  "
         //           << value << std::endl;
      }
      assert(paramValues[i] == value);
      assert(m_func->getParamValue(paramNames[i]) == value);
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

void FunctionTest::freeParameters(const std::vector<Parameter> &params) {

// Create a local copy of the test Parameters
   std::vector<Parameter> my_params = params;

// Make every other parameter free, the others fixed
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if ((i % 2) == 0) {
         my_params[i].setFree(true);
      } else {
         my_params[i].setFree(false);
      }
   }

// Sync up the parameter names
   std::vector<std::string> paramNames;
   m_func->getParamNames(paramNames);
   my_params.resize(paramNames.size());
   for (unsigned int i = 0; i < my_params.size(); i++) {
      my_params[i].setName(paramNames[i]);
   }

// set the all the Parameters
   m_func->setParams(my_params);

// retrieve the free ones
   std::vector<std::string> freeParamNames;
   m_func->getFreeParamNames(freeParamNames);
   std::vector<double> freeParamValues;
   m_func->getFreeParamValues(freeParamValues);

// and check against my_params
   int j = 0;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (my_params[i].isFree()) {
         assert(my_params[i].getName() == freeParamNames[j]);
         assert(my_params[i].getValue() == freeParamValues[j]);
         j++;
      }
   }

// reset the free Parameter values as a group
   for (unsigned int i = 0; i < freeParamValues.size(); i++) {
      freeParamValues[i] /= 5.;
   }
   m_func->setFreeParamValues(freeParamValues);

// check
   j = 0;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (my_params[i].isFree()) {
         assert(my_params[i].getValue()/5. 
                == m_func->getParamValue(freeParamNames[j]));
         j++;
      }
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

void FunctionTest::funcEvaluations(const std::vector<Arg*> &arguments,
                                   const std::vector<double> &returnValues) {
   if (arguments.size() != returnValues.size()) {
      throw Exception(
         std::string("FunctionTest::funcEvaluations:\n") + 
         std::string("Sizes of arguments and (expected) returnValues") + 
         std::string(" do not match."));
   }

//    std::vector<double> params;
//    m_func->getParamValues(params);
//    for (unsigned int i = 0; i < params.size(); i++) 
//       std::cout << params[i] << std::endl;
   for (unsigned int i = 0; i < arguments.size(); i++) {
      double val = m_func->value(*arguments[i]);
      assert(val == returnValues[i]);
   }
}

void FunctionTest::derivatives(const std::vector<Arg*> &arguments,
                               double eps) {
// loop over all Parameters
   for (unsigned int k = 0; k < arguments.size(); k++) {
      Arg *my_arg = arguments[k];

// Check the derivatives wrt all Parameters against numerical estimates 
      std::vector<double> derivs;
      m_func->getFreeDerivs(*my_arg, derivs);
      std::vector<double> params;
      m_func->getFreeParamValues(params);
      for (size_t i(0); i < params.size(); i++) {
         // std::cout << i << "  "
         //           << params[i] << std::endl;
      }

      double f0 = m_func->value(*my_arg);
      // std::cout << dynamic_cast<optimizers::dArg *>(my_arg)->getValue() << "  " 
      //           << f0 << std::endl;
      for (unsigned int i = 0; i < params.size(); i++) {
         std::vector<double> new_params = params;
         double delta1 = new_params[i]*eps;
         delta1 = ::select_delta(new_params[i], delta1);
         new_params[i] += delta1;
         m_func->setFreeParamValues(new_params);
         double f1 = m_func->value(*my_arg);

         new_params = params;
         double delta2 = -new_params[i]*eps;
         delta2 = ::select_delta(new_params[i], delta2);
         new_params[i] += delta2;
         m_func->setFreeParamValues(new_params);
         double f2 = m_func->value(*my_arg);

         double my_deriv = (f1 - f2)/(delta1 - delta2);

//          std::cout << i << "  "
//                    << params[i] << "  "
//                    << f1 << "  "
//                    << f2 << "  "
//                    << (f1 - f2) << "  "
//                    << delta1 << "  "
//                    << delta2 << "  ";

//          std::cout << my_deriv << "  "
//                    << derivs[i] << "  ";
         if (derivs[i] != 0) {
            double value = fabs(my_deriv/derivs[i] - 1.);
//             std::cout << value << std::endl;
            assert(value < eps*10.);
         } else {
//             std::cout << 0 << std::endl;
            assert(fabs(my_deriv) < eps*10.);
         }
      }

// Check that the free derivatives are being accessed correctly.
      m_func->setFreeParamValues(params);     // restore the Parameter values
      std::vector<Parameter> parameters;  // get all Parameter info
      m_func->getParams(parameters);
      std::vector<double> freeDerivs;     // and the free derivatives
      m_func->getFreeDerivs(*my_arg, freeDerivs);
      int j = 0;
      for (unsigned int i = 0; i < parameters.size(); i++) {
         if (parameters[i].isFree()) {
            // The var1, var2 and the fpe test are needed to ensure
            // that these variables are accessed properly in the
            // assert for rh9_gcc32opt and rhel4_gcc34opt builds.
            double var1(m_func->derivByParam(*my_arg, parameters[i].getName()));
            double var2(var1*freeDerivs[j]);
            if (var2 != var2) {
               throw std::runtime_error("FPE in FunctionTest::derivatives");
            }
            double x1 = m_func->derivByParam(*my_arg, parameters[i].getName());
	    double x2 = freeDerivs[j];
	    if (x1 != 0) {
	      double dx = fabs((x1 - x2) / x1);
	      assert(dx <= 1e-6);
	    } else if (x2 != 0) {
	      double dx = fabs((x1 - x2) / x2);
	      assert(dx <= 1e-6);
	    } else {
	      // x1 == x2 == 0
	    }
            j++;
         }
      }
   }

// restore the original Parameters
   m_func->setParams(m_originalParameters);
}

} //namespace optimizers
