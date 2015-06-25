/**
 * @file FunctionFactory.h
 * @brief Prototype pattern for generating Function objects.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/FunctionFactory.h,v 1.5 2015/03/19 18:31:25 jchiang Exp $
 */

#ifndef optimizers_FunctionFactory_h
#define optimizers_FunctionFactory_h

#include <map>
#include <string>

#include "optimizers/Function.h"
#include "optimizers/Exception.h"

namespace optimizers {

/**
 * @class FunctionFactory
 *
 * @brief Use the prototype pattern to supply customized Function
 * objects, the parameters of which are specified by an xml input
 * file.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/FunctionFactory.h,v 1.5 2015/03/19 18:31:25 jchiang Exp $
 */

class FunctionFactory {

public:

   FunctionFactory();

   virtual ~FunctionFactory();

   Function * create(const std::string & name);

   void addFunc(const std::string & name, Function * func,
                bool fromClone=true);

   void addFunc(Function * func, bool fromClone=true);

   void getFunctionNames(std::vector<std::string> & funcNames);

   void readXml(const std::string & xmlFile);

   void writeXml(const std::string & outputFile);

private:

   std::map<std::string, Function *> m_prototypes;

};

} // namespace optimizers

#endif // optimizers_FunctionFactory_h
