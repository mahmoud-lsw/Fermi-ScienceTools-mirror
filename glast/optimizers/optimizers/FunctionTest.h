/**
 * @file FunctionTest.h
 * @brief Declaration of unit test code for Function class hierarchy
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/FunctionTest.h,v 1.3 2004/12/30 00:24:40 jchiang Exp $
 */


#ifndef optimizers_FunctionTest_h
#define optimizers_FunctionTest_h

#include <cassert>

#include "optimizers/Parameter.h"
#include "optimizers/Arg.h"
#include "optimizers/Function.h"
#include "optimizers/Exception.h"

namespace optimizers {

/**
 * @class FunctionTest
 * @brief Unit test code for Function class hierarchy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/FunctionTest.h,v 1.3 2004/12/30 00:24:40 jchiang Exp $
 */

class FunctionTest {

public:

   FunctionTest(Function &func, const std::string &name) : m_func(&func) {
      m_func->setName(name);
      assert(m_func->getName() == name);
      m_func->getParams(m_originalParameters);
   }

   ~FunctionTest() {}

   void parameters(const std::vector<Parameter> & params);

   void freeParameters(const std::vector<Parameter> & params);

   void funcEvaluations(const std::vector<Arg *> & arguments,
                        const std::vector<double> & returnValues);

   void derivatives(const std::vector<Arg *> & arguments,
                    double eps = 1e-5);
      
private:

   Function * m_func;
   std::vector<Parameter> m_originalParameters;

};

} // namespace FunctionTest

#endif // optimizers_FunctionTest_h
