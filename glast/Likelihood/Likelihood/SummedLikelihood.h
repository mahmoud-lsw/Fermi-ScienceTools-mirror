/**
 * @file SummedLikelihood.h
 *
 * @brief Statistic object that comprises LogLike objects that have
 * identical source models.
 * 
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/SummedLikelihood.h,v 1.1.1.3.6.5 2015/04/26 06:53:44 jasercio Exp $
 */

#ifndef Likelihood_SummedLikelihood_h
#define Likelihood_SummedLikelihood_h

#include <map>
#include <set>
#include <string>
#include <vector>

#include "optimizers/Statistic.h"

#include "Likelihood/LogLike.h"
#include "Likelihood/TiedParameter.h"

namespace Likelihood {

/*
 * @class SummedLikelihood
 */

class SummedLikelihood : public optimizers::Statistic {

public:

   SummedLikelihood() : optimizers::Statistic(), m_masterComponent(0) {}

   virtual ~SummedLikelihood() throw();

   void addComponent(LogLike & component);

   virtual double value() const;
   virtual void getFreeParams(std::vector<optimizers::Parameter> &params) const;
   virtual void setFreeParamValues(const std::vector<double> & values);
   virtual unsigned int getNumFreeParams() const;
   virtual void getFreeDerivs(std::vector<double> & derivs) const;

   void syncParams();

   double NpredValue(const std::string &) const {return 0;}

   /// Member functions to support tying of parameters.

   void tieParameters(const std::vector<size_t> & par_indices);
   void setErrors(const std::vector<double> & errors);

   /// @return The index used by Minos for free parameters.
   int findIndex(size_t par_index) const;

   /// @return The TiedParameter object containing par_index.
   /// @param par_index The index used by SourceModel for both free
   ///        and fixed parameters.
   TiedParameter & getTiedParam(size_t par_index);
   
   void setTiedParamValue(size_t par_index, double value);

protected:

   double value(optimizers::Arg &) const {
      return value();
   }

   double derivByParamImp(optimizers::Arg &, const std::string &) const {
      return 0;
   }

   optimizers::Function * clone() const {return 0;}

   virtual void fetchParamValues(std::vector<double> & values,
                                 bool getFree) const;

private:

   typedef std::vector<LogLike *> ComponentVector_t;
   typedef ComponentVector_t::iterator ComponentIterator_t;
   typedef ComponentVector_t::const_iterator ComponentConstIterator_t;
   ComponentVector_t m_components;

   std::vector<TiedParameter *> m_tiedPars;

   /// Set of all parameter indices that are associated with tied parameters.
   std::set<size_t> m_tiedIndices;

   LogLike * m_masterComponent;

};

} // namespace Likelihood

#endif // Likelihood_SummedLikelihood_h