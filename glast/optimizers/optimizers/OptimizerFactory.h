/**
 * @file OptimizerFactory.h
 * @brief Factory for generating Optimizer objects by name.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/OptimizerFactory.h,v 1.1 2006/05/10 18:06:26 jchiang Exp $
 */

#ifndef optimizer_OptimizerFactory_h
#define optimizer_OptimizerFactory_h

#include <string>

namespace optimizers {

class Optimizer;
class Statistic;

/**
 * @class OptimizerFactory
 * @brief Singleton factory for generating Optimizer objects by name.
 */

class OptimizerFactory {

public:

   static OptimizerFactory & instance();

   Optimizer * create(const std::string & optimizerName,
                      Statistic & stat);

protected:

   OptimizerFactory() {}

private:

   static OptimizerFactory * s_instance;

};

} // namespace optimizers

#endif // optimizer_OptimizerFactory_h
