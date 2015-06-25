/**
 * @file OptimizerFactory.h
 * @brief Implementation for Singleton factory to create Optimizer
 * objects by name.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/OptimizerFactory.cxx,v 1.8 2009/05/27 20:53:10 pln Exp $
 */

#include <stdexcept>

#include "optimizers/Drmngb.h"
#include "optimizers/Drmnfb.h"
#include "optimizers/ModNewton.h"
#include "optimizers/Lbfgs.h"
#include "optimizers/Minuit.h"
#include "optimizers/NewMinuit.h"
#include "optimizers/Optimizer.h"
#include "optimizers/Statistic.h"
#include "optimizers/Powell.h"

#include "optimizers/OptimizerFactory.h"

namespace optimizers {

OptimizerFactory * OptimizerFactory::s_instance(0);

OptimizerFactory & OptimizerFactory::instance() {
   if (s_instance == 0) {
      s_instance = new OptimizerFactory();
   }
   return *s_instance;
}

Optimizer * OptimizerFactory::create(const std::string & optimizerName,
                                     Statistic & stat) {
   if (optimizerName == "Minuit" || optimizerName == "MINUIT") {
      return new Minuit(stat);
   } else if (optimizerName == "Lbfgs" || optimizerName == "LBFGS") {
      return new Lbfgs(stat);
   } else if (optimizerName == "Drmngb" || optimizerName == "DRMNGB") {
//      return new Drmngb(stat);
      return new ModNewton(stat, true);
   } else if (optimizerName == "Drmnfb" || optimizerName == "DRMNFB") {
//      return new Drmnfb(stat);
      return new ModNewton(stat, false);
   } else if (optimizerName == "Powell" || optimizerName == "POWELL") {
      return new Powell(stat);
   } else if (optimizerName == "NewMinuit" || optimizerName == "NEWMINUIT") {
      return new NewMinuit(stat);
   } else {
      throw std::runtime_error("Invalid optimizer choice: " + optimizerName);
   }
}

} // namespace optimizers
