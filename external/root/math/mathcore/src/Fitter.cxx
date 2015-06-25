// @(#)root/mathcore:$Id: Fitter.cxx,v 1.1.1.3 2013/01/23 16:01:36 areustle Exp $
// Author: L. Moneta Mon Sep  4 17:00:10 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class Fitter


#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/FcnAdapter.h"
#include "Math/Error.h"

#include <memory> 

#include "Math/IParamFunction.h" 

#include "Math/MultiDimParamFunctionAdapter.h"

// #include "TMatrixDSym.h"
// for debugging
//#include "TMatrixD.h"
// #include <iomanip>

namespace ROOT { 

   namespace Fit { 

// use a static variable to get default minimizer options for error def
// to see if user has changed it later on. If it has not been changed we set 
// for the likelihood method an error def of 0.5 
// t.b.d : multiply likelihood by 2 so have same error def definition as chi2 
      double gDefaultErrorDef = ROOT::Math::MinimizerOptions::DefaultErrorDef();


Fitter::Fitter() : 
   fUseGradient(false),
   fBinFit(false),
   fFitType(0),
   fDataSize(0),
   fFunc(0)
{
   // Default constructor implementation.
   fResult = std::auto_ptr<ROOT::Fit::FitResult>(new ROOT::Fit::FitResult() );
}

Fitter::~Fitter() 
{
   // Destructor implementation.
   // delete function if not empty
   if (fFunc) delete fFunc; 
}

Fitter::Fitter(const Fitter & rhs) 
{
   // Implementation of copy constructor.
   // copy FitResult, FitConfig and clone fit function
   (*this) = rhs; 
}

Fitter & Fitter::operator = (const Fitter &rhs) 
{
   // Implementation of assignment operator.
   // dummy implementation, since it is private
   if (this == &rhs) return *this;  // time saving self-test
//    fUseGradient = rhs.fUseGradient; 
//    fBinFit = rhs.fBinFit; 
//    fResult = rhs.fResult;
//    fConfig = rhs.fConfig; 
//    // function is copied and managed by FitResult (maybe should use an auto_ptr)
//    fFunc = fResult.ModelFunction(); 
//    if (rhs.fFunc != 0 && fResult.ModelFunction() == 0) { // case no fit has been done yet - then clone 
//       if (fFunc) delete fFunc; 
//       fFunc = dynamic_cast<IModelFunction *>( (rhs.fFunc)->Clone() ); 
//       assert(fFunc != 0); 
//    }
   return *this; 
}

void Fitter::SetFunction(const IModelFunction & func, bool useGradient) 
{

   fUseGradient = useGradient; 
   if (fUseGradient) { 
      const IGradModelFunction * gradFunc = dynamic_cast<const IGradModelFunction*>(&func);
      if (gradFunc) { 
         SetFunction(*gradFunc, true);
         return;
      }
      else {
         MATH_WARN_MSG("Fitter::SetFunction","Requested function does not provide gradient - use it as non-gradient function ");
      }
   }
   fUseGradient = false;
 
   //  set the fit model function (clone the given one and keep a copy ) 
   //std::cout << "set a non-grad function" << std::endl; 

   fFunc = dynamic_cast<IModelFunction *>(func.Clone() ); 
   assert(fFunc != 0);
   
   // creates the parameter  settings 
   fConfig.CreateParamsSettings(*fFunc); 
}


void Fitter::SetFunction(const IModel1DFunction & func, bool useGradient) 
{ 
   fUseGradient = useGradient; 
   if (fUseGradient) { 
      const IGradModel1DFunction * gradFunc = dynamic_cast<const IGradModel1DFunction*>(&func);
      if (gradFunc) { 
         SetFunction(*gradFunc, true);
         return;
      }
      else {
         MATH_WARN_MSG("Fitter::SetFunction","Requested function does not provide gradient - use it as non-gradient function ");
      }
   }
   fUseGradient = false;
   //std::cout << "set a 1d function" << std::endl; 

   // function is cloned when creating the adapter
   fFunc = new ROOT::Math::MultiDimParamFunctionAdapter(func);

   // creates the parameter  settings 
   fConfig.CreateParamsSettings(*fFunc); 
}

void Fitter::SetFunction(const IGradModelFunction & func, bool useGradient) 
{ 
   fUseGradient = useGradient;
   //std::cout << "set a grad function" << std::endl; 
   //  set the fit model function (clone the given one and keep a copy ) 
   fFunc = dynamic_cast<IGradModelFunction *> ( func.Clone() ); 
   assert(fFunc != 0);

   // creates the parameter  settings 
   fConfig.CreateParamsSettings(*fFunc); 
}


void Fitter::SetFunction(const IGradModel1DFunction & func, bool useGradient) 
{ 
   //std::cout << "set a 1d grad function" << std::endl; 
   fUseGradient = useGradient;
   // function is cloned when creating the adapter
   fFunc = new ROOT::Math::MultiDimParamGradFunctionAdapter(func);

   // creates the parameter  settings 
   fConfig.CreateParamsSettings(*fFunc); 
}


bool Fitter::SetFCN(const ROOT::Math::IMultiGenFunction & fcn, const double * params, unsigned int dataSize, bool chi2fit) {
   // set the objective function for the fit 
   // if params is not NULL create the parameter settings
   fUseGradient = false; 
   unsigned int npar  = fcn.NDim(); 
   if (npar == 0) { 
      MATH_ERROR_MSG("Fitter::SetFCN","FCN function has zero parameters ");
      return false;
   }
   if (params != 0 ) 
      fConfig.SetParamsSettings(npar, params);
   else {
      if ( fConfig.ParamsSettings().size() != npar) { 
         MATH_ERROR_MSG("Fitter::SetFCN","wrong fit parameter settings");
         return false;
      }
   }

   fBinFit = chi2fit; 
   fDataSize = dataSize; 

   // keep also a copy of FCN function and set this in minimizer so they will be managed together
   // (remember that cloned copy will still depends on data and model function pointers) 
   fObjFunction = std::auto_ptr<ROOT::Math::IMultiGenFunction> ( fcn.Clone() ); 

   return true;
}

bool Fitter::SetFCN(const ROOT::Math::IMultiGradFunction & fcn, const double * params, unsigned int dataSize , bool chi2fit) {
   // set the objective function for the fit 
   // if params is not NULL create the parameter settings
   if (!SetFCN(static_cast<const ROOT::Math::IMultiGenFunction &>(fcn),params, dataSize, chi2fit) ) return false;
   fUseGradient = true; 
   return true;
}

bool Fitter::SetFCN(const ROOT::Math::FitMethodFunction & fcn, const double * params) {
   // set the objective function for the fit 
   // if params is not NULL create the parameter settings
   bool chi2fit = (fcn.Type() == ROOT::Math::FitMethodFunction::kLeastSquare );
   if (!SetFCN(fcn,params,fcn.NPoints(), chi2fit) ) return false;
   fUseGradient = false; 
   fFitType = fcn.Type();
   return true;
}

bool Fitter::SetFCN(const ROOT::Math::FitMethodGradFunction & fcn, const double * params) {
   // set the objective function for the fit 
   // if params is not NULL create the parameter settings
   bool chi2fit = (fcn.Type() == ROOT::Math::FitMethodGradFunction::kLeastSquare );
   if (!SetFCN(fcn,params,fcn.NPoints(), chi2fit) ) return false;
   fUseGradient = true; 
   fFitType = fcn.Type();
   return true;
}



bool Fitter::FitFCN(const BaseFunc & fcn, const double * params, unsigned int dataSize, bool chi2fit) { 
   // fit a user provided FCN function
   // create fit parameter settings
   if (!SetFCN(fcn, params,dataSize,chi2fit) ) return false;
   return FitFCN(); 
}



bool Fitter::FitFCN(const BaseGradFunc & fcn, const double * params, unsigned int dataSize, bool chi2fit) { 
   // fit a user provided FCN gradient function

   if (!SetFCN(fcn, params,dataSize, chi2fit) ) return false;
   return FitFCN();
}

bool Fitter::FitFCN(const ROOT::Math::FitMethodFunction & fcn, const double * params) {
   // fit using the passed objective function for the fit 
   if (!SetFCN(fcn, params) ) return false;
   return FitFCN();
}

bool Fitter::FitFCN(const ROOT::Math::FitMethodGradFunction & fcn, const double * params) {
   // fit using the passed objective function for the fit 
   if (!SetFCN(fcn, params) ) return false;
   return FitFCN();
}


bool Fitter::SetFCN(MinuitFCN_t fcn, int npar, const double * params , unsigned int dataSize , bool chi2fit ){ 
   // set TMinuit style FCN type (global function pointer)  
   // create corresponfing objective function from that function

   if (npar == 0) {
      npar = fConfig.ParamsSettings().size(); 
      if (npar == 0) { 
         MATH_ERROR_MSG("Fitter::FitFCN","Fit Parameter settings have not been created ");
         return false;
      }
   }

   ROOT::Fit::FcnAdapter  newFcn(fcn,npar); 
   return SetFCN(newFcn,params,dataSize,chi2fit); 
}

bool Fitter::FitFCN(MinuitFCN_t fcn, int npar, const double * params , unsigned int dataSize , bool chi2fit ) { 
   // fit using Minuit style FCN type (global function pointer)  
   // create corresponfing objective function from that function
   if (!SetFCN(fcn, npar, params, dataSize, chi2fit)) return false; 
   fUseGradient = false; 
   return FitFCN(); 
}

bool Fitter::FitFCN() { 
   // fit using the previously set  FCN function

   // in case a model function exists from a previous fit delete it
   if (fFunc && fResult->FittedFunction() == 0) delete fFunc; 
   fFunc = 0;

   if (fObjFunction.get() == 0) { 
      MATH_ERROR_MSG("Fitter::FitFCN","Objective function has not been set");
      return false; 
   }
   // init the minimizer
   if (!DoInitMinimizer() ) return false; 
   // perform the minimization
   return DoMinimization(); 
}

bool Fitter::EvalFCN() { 
   // evaluate the FCN using the stored values in fConfig

   if (fFunc && fResult->FittedFunction() == 0) delete fFunc; 
   fFunc = 0;

   if (fObjFunction.get() == 0) { 
      MATH_ERROR_MSG("Fitter::FitFCN","Objective function has not been set");
      return false; 
   }
   // create a Fit result from the fit configuration ); 
   fResult = std::auto_ptr<ROOT::Fit::FitResult>(new ROOT::Fit::FitResult(fConfig) );
   // evaluate one time the FCN 
   double fcnval = (*fObjFunction)(fResult->GetParams() );
   // update fit result
   fResult->fVal = fcnval;
   fResult->fNCalls++;
   return true;
}


bool Fitter::DoLeastSquareFit(const BinData & data) { 
   // perform a chi2 fit on a set of binned data 
   
   // check function
   if (fFunc == 0) { 
      MATH_ERROR_MSG("Fitter::DoLeastSquareFit","model function is not set");
      return false; 
   }

#ifdef DEBUG
   std::cout << "Fitter ParamSettings " << Config().ParamsSettings()[3].IsBound() << " lower limit " <<  Config().ParamsSettings()[3].LowerLimit() << " upper limit " <<  Config().ParamsSettings()[3].UpperLimit() << std::endl;
#endif

   fBinFit = true; 
   fDataSize = data.Size();

   // check if fFunc provides gradient
   if (!fUseGradient) { 
      // do minimzation without using the gradient
      Chi2FCN<BaseFunc> chi2(data,*fFunc); 
      fFitType = chi2.Type();
      return DoMinimization (chi2); 
   } 
   else { 
      // use gradient 
      if (fConfig.MinimizerOptions().PrintLevel() > 0) 
         MATH_INFO_MSG("Fitter::DoLeastSquareFit","use gradient from model function");        
      IGradModelFunction * gradFun = dynamic_cast<IGradModelFunction *>(fFunc); 
      if (gradFun != 0) { 
         Chi2FCN<BaseGradFunc> chi2(data,*gradFun); 
         fFitType = chi2.Type();
         return DoMinimization (chi2); 
      }
      MATH_ERROR_MSG("Fitter::DoLeastSquareFit","wrong type of function - it does not provide gradient");
   }
   return false; 
}

bool Fitter::DoLikelihoodFit(const BinData & data, bool extended) { 
   // perform a likelihood fit on a set of binned data 
   // The fit is extended (Poisson logl_ by default 
   
   bool  useWeight = fConfig.UseWeightCorrection();  

   // check function
   if (fFunc == 0) { 
      MATH_ERROR_MSG("Fitter::DoLikelihoodFit","model function is not set");
      return false; 
   }

   // logl fit (error should be 0.5) set if different than default values (of 1)
   if (fConfig.MinimizerOptions().ErrorDef() == gDefaultErrorDef ) { 
      fConfig.MinimizerOptions().SetErrorDef(0.5);
   }

   if (useWeight && fConfig.MinosErrors() ) {
      MATH_INFO_MSG("Fitter::DoLikelihoodFit","MINOS errors cannot be computed in weighted likelihood fits");
      fConfig.SetMinosErrors(false);
   }
      

   fBinFit = true; 
   fDataSize = data.Size();

   // create a chi2 function to be used for the equivalent chi-square
   Chi2FCN<BaseFunc> chi2(data,*fFunc); 

   if (!fUseGradient) { 
      // do minimization without using the gradient
      PoissonLikelihoodFCN<BaseFunc> logl(data,*fFunc, useWeight, extended); 
      fFitType = logl.Type();
      // do minimization
      if (!DoMinimization (logl, &chi2) ) return false; 
      if (useWeight) { 
         logl.UseSumOfWeightSquare();
         if (!ApplyWeightCorrection(logl) ) return false; 
      }
   } 
   else { 
      if (fConfig.MinimizerOptions().PrintLevel() > 0) 
         MATH_INFO_MSG("Fitter::DoLikelihoodFit","use gradient from model function");        
      // check if fFunc provides gradient
      IGradModelFunction * gradFun = dynamic_cast<IGradModelFunction *>(fFunc); 
      if (gradFun == 0) { 
         MATH_ERROR_MSG("Fitter::DoLikelihoodFit","wrong type of function - it does not provide gradient");        
         return false; 
      }
      // use gradient for minimization
      // not-extended is not impelemented in this case 
      if (!extended) {  
         MATH_WARN_MSG("Fitter::DoLikelihoodFit","Not-extended binned fit with gradient not yet supported - do an extended fit");        
      }
      PoissonLikelihoodFCN<BaseGradFunc> logl(data,*gradFun, useWeight, true); 
      fFitType = logl.Type();
      // do minimization
      if (!DoMinimization (logl, &chi2) ) return false;
      if (useWeight) { 
         logl.UseSumOfWeightSquare();
         if (!ApplyWeightCorrection(logl) ) return false; 
      }
   }
   return true; 
}


bool Fitter::DoLikelihoodFit(const UnBinData & data, bool extended) { 
   // perform a likelihood fit on a set of unbinned data 

   bool useWeight = fConfig.UseWeightCorrection(); 

   if (fFunc == 0) { 
      MATH_ERROR_MSG("Fitter::DoLikelihoodFit","model function is not set");
      return false; 
   }

   if (useWeight && fConfig.MinosErrors() ) {
      MATH_INFO_MSG("Fitter::DoLikelihoodFit","MINOS errors cannot be computed in weighted likelihood fits");
      fConfig.SetMinosErrors(false);
   }


   fBinFit = false; 
   fDataSize = data.Size();

#ifdef DEBUG
   int ipar = 0;
   std::cout << "Fitter ParamSettings " << Config().ParamsSettings()[ipar].IsBound() << " lower limit " <<  Config().ParamsSettings()[ipar].LowerLimit() << " upper limit " <<  Config().ParamsSettings()[ipar].UpperLimit() << std::endl;
#endif

   // logl fit (error should be 0.5) set if different than default values (of 1)
   if (fConfig.MinimizerOptions().ErrorDef() == gDefaultErrorDef ) {
      fConfig.MinimizerOptions().SetErrorDef(0.5);
   }

   if (!fUseGradient) { 
      // do minimization without using the gradient
      LogLikelihoodFCN<BaseFunc> logl(data,*fFunc, useWeight, extended); 
      fFitType = logl.Type();
      if (!DoMinimization (logl) ) return false;
      if (useWeight) { 
         logl.UseSumOfWeightSquare();
         if (!ApplyWeightCorrection(logl) ) return false; 
      }
      return true;
   } 
   else { 
      // use gradient : check if fFunc provides gradient
      if (fConfig.MinimizerOptions().PrintLevel() > 0) 
         MATH_INFO_MSG("Fitter::DoLikelihoodFit","use gradient from model function");        
      IGradModelFunction * gradFun = dynamic_cast<IGradModelFunction *>(fFunc); 
      if (gradFun != 0) { 
         if (extended) {  
            MATH_WARN_MSG("Fitter::DoLikelihoodFit","Extended unbinned fit with gradient not yet supported - do a not-extended fit");        
         }
         LogLikelihoodFCN<BaseGradFunc> logl(data,*gradFun,useWeight, extended); 
         fFitType = logl.Type();
         if (!DoMinimization (logl) ) return false;
         if (useWeight) { 
            logl.UseSumOfWeightSquare();
            if (!ApplyWeightCorrection(logl) ) return false; 
         }
         return true;
      }
      MATH_ERROR_MSG("Fitter::DoLikelihoodFit","wrong type of function - it does not provide gradient");
   }      
   return false; 
}


bool Fitter::DoLinearFit(const BinData & data ) { 

   // perform a linear fit on a set of binned data 
   std::string  prevminimizer = fConfig.MinimizerType();  
   fConfig.SetMinimizer("Linear"); 

   fBinFit = true; 

   bool ret =  DoLeastSquareFit(data); 
   fConfig.SetMinimizer(prevminimizer.c_str());
   return ret; 
}


bool Fitter::CalculateHessErrors() { 
   // compute the Hesse errors according to configuration
   // set in the parameters and append value in fit result
   if (fObjFunction.get() == 0) { 
      MATH_ERROR_MSG("Fitter::CalculateHessErrors","Objective function has not been set");
      return false; 
   }
   // assume a fResult pointer always exists
   assert (fResult.get() );

   // need a special treatment in case of weighted likelihood fit 
   // (not yet implemented) 
   if (fFitType == 2 && fConfig.UseWeightCorrection() ) { 
      MATH_ERROR_MSG("Fitter::CalculateHessErrors","Re-computation of Hesse errors not implemented for weighted likelihood fits");
      MATH_INFO_MSG("Fitter::CalculateHessErrors","Do the Fit using configure option FitConfig::SetParabErrors()");
      return false;
   }
      // if (!fUseGradient ) {
      // ROOT::Math::FitMethodFunction * fcn = dynamic_cast< ROOT::Math::FitMethodFunction *>(fObjFunction.get());
      // if (fcn  && fcn->Type() ==  ROOT::Math::FitMethodFunction::kLogLikelihood) { 
      //    if (!fBinFit) { 
      //       ROOT::Math::LogLikelihoodFunction * nll = dynamic_cast< ROOT::Math::LogLikelihoodFunction *>(fcn);
      //       assert(nll); 
      //       nll->UseSumOfWeightSquare(false);
      //    } 
      //    else { 
      //       ROOT::Math::PoissonLikelihoodFunction * nll = dynamic_cast< ROOT::Math::PoissonLikelihoodFunction *>(fcn);
      //       assert(nll); 
      //       nll->UseSumOfWeightSquare(false);
      //    }   
      //    // reset fcn in minimizer
      // }


   // create minimizer if not done or if name has changed 
   if ( !fMinimizer.get()  || 
       fResult->MinimizerType().find(fConfig.MinimizerType()) == std::string::npos ) {
      bool ret = DoInitMinimizer(); 
      if (!ret) {
       MATH_ERROR_MSG("Fitter::CalculateHessErrors","Error initializing the minimizer");
       return false; 
      }
   }
      

   if (!fMinimizer.get() ) { 
       MATH_ERROR_MSG("Fitter::CalculateHessErrors","Need to do a fit before calculating the errors");
       return false; 
   }

   //run Hesse
   bool ret = fMinimizer->Hesse();
   if (!ret) MATH_WARN_MSG("Fitter::CalculateHessErrors","Error when calculating Hessian");


   // update minimizer results with what comes out from Hesse
   // in case is empty - create from a FitConfig
   if (fResult->IsEmpty() )    
      fResult = std::auto_ptr<ROOT::Fit::FitResult>(new ROOT::Fit::FitResult(fConfig) );


   ret |= fResult->Update(*fMinimizer, ret); 

   // when possible get ncalls from FCN and set in fit result
   if (fFitType != ROOT::Math::FitMethodFunction::kUndefined ) { 
      fResult->fNCalls = GetNCallsFromFCN(); 
   }

   // set also new errors in FitConfig
   if (fConfig.UpdateAfterFit() && ret) DoUpdateFitConfig(); 

   return ret; 
}


bool Fitter::CalculateMinosErrors() { 
   // compute the Minos errors according to configuration
   // set in the parameters and append value in fit result
   
   // in case it has not been set - set by default Minos on all parameters 
   fConfig.SetMinosErrors(); 

   if (!fMinimizer.get() ) { 
       MATH_ERROR_MSG("Fitter::CalculateMinosErrors","Minimizer does not exist - cannot calculate Minos errors");
       return false; 
   }

   if (!fResult.get() || fResult->IsEmpty() ) { 
       MATH_ERROR_MSG("Fitter::CalculateMinosErrors","Invalid Fit Result - cannot calculate Minos errors");
       return false; 
   }

   if (fFitType == 2 && fConfig.UseWeightCorrection() ) { 
      MATH_ERROR_MSG("Fitter::CalculateMinosErrors","Computation of MINOS errors not implemented for weighted likelihood fits");
      return false;
   }


   const std::vector<unsigned int> & ipars = fConfig.MinosParams(); 
   unsigned int n = (ipars.size() > 0) ? ipars.size() : fResult->Parameters().size(); 
   bool ok = false; 
   for (unsigned int i = 0; i < n; ++i) {
      double elow, eup;
      unsigned int index = (ipars.size() > 0) ? ipars[i] : i; 
      bool ret = fMinimizer->GetMinosError(index, elow, eup);
      if (ret) fResult->SetMinosError(index, elow, eup); 
      ok |= ret; 
   }
   if (!ok) 
       MATH_ERROR_MSG("Fitter::CalculateMinosErrors","Minos error calculation failed for all parameters");

   return ok; 
}



// traits for distinhuishing fit methods functions from generic objective functions
template<class Func> 
struct ObjFuncTrait { 
   static unsigned int NCalls(const Func &  ) { return 0; }
   static int Type(const Func & ) { return -1; }
   static bool IsGrad() { return false; }
};
template<>
struct ObjFuncTrait<ROOT::Math::FitMethodFunction> { 
   static unsigned int NCalls(const ROOT::Math::FitMethodFunction & f ) { return f.NCalls(); }
   static int Type(const ROOT::Math::FitMethodFunction & f) { return f.Type(); }
   static bool IsGrad() { return false; }
};
template<>
struct ObjFuncTrait<ROOT::Math::FitMethodGradFunction> { 
   static unsigned int NCalls(const ROOT::Math::FitMethodGradFunction & f ) { return f.NCalls(); }
   static int Type(const ROOT::Math::FitMethodGradFunction & f) { return f.Type(); }
   static bool IsGrad() { return true; }
};

bool Fitter::DoInitMinimizer() { 
   //initialize minimizer by creating it
   // and set there the objective function 
   // obj function must have been copied before 
   assert(fObjFunction.get() );

   // check configuration and objective  function
   if ( fConfig.ParamsSettings().size() != fObjFunction->NDim() ) { 
      MATH_ERROR_MSG("Fitter::DoInitMinimizer","wrong function dimension or wrong size for FitConfig");
      return false;
   }

   // create first Minimizer  
   // using an auto_Ptr will delete the previous existing one
   fMinimizer = std::auto_ptr<ROOT::Math::Minimizer> ( fConfig.CreateMinimizer() );
   if (fMinimizer.get() == 0) {
      MATH_ERROR_MSG("Fitter::FitFCN","Minimizer cannot be created");
      return false; 
   }

   // in case of gradient function one needs to downcast the pointer
   if (fUseGradient) { 
      const ROOT::Math::IMultiGradFunction * gradfcn = dynamic_cast<const ROOT::Math::IMultiGradFunction *> (fObjFunction.get() );
      if (!gradfcn) { 
         MATH_ERROR_MSG("Fitter::DoInitMinimizer","wrong type of function - it does not provide gradient");
         return false;
      } 
      fMinimizer->SetFunction( *gradfcn);
   }
   else 
      fMinimizer->SetFunction( *fObjFunction);


   fMinimizer->SetVariables(fConfig.ParamsSettings().begin(), fConfig.ParamsSettings().end() ); 

   // if requested parabolic error do correct error  analysis by the minimizer (call HESSE) 
   if (fConfig.ParabErrors()) fMinimizer->SetValidError(true);

   return true; 

}


bool Fitter::DoMinimization(const ROOT::Math::IMultiGenFunction * chi2func) { 
   // perform the minimization (assume we have already initialized the minimizer)

   assert(fMinimizer.get() );
   
   bool ret = fMinimizer->Minimize(); 

   // unsigned int ncalls =  ObjFuncTrait<ObjFunc>::NCalls(*fcn);
   // int fitType =  ObjFuncTrait<ObjFunc>::Type(objFunc);


   fResult = std::auto_ptr<FitResult> ( new FitResult(*fMinimizer,fConfig, fFunc, ret, fDataSize, 
                                                      fBinFit, chi2func ) );

   // when possible get ncalls from FCN and set in fit result
   if (fResult->fNCalls == 0 && fFitType != ROOT::Math::FitMethodFunction::kUndefined ) { 
      fResult->fNCalls = GetNCallsFromFCN(); 
   }


#ifdef DEBUG
      std::cout << "ROOT::Fit::Fitter::DoMinimization : ncalls = " << fResult->fNCalls << " type of objfunc " << fFitFitResType << "  typeid: " << typeid(*fObjFunction).name() << " use gradient " << fUseGradient << std::endl;
#endif

   if (fConfig.NormalizeErrors() && fFitType == ROOT::Math::FitMethodFunction::kLeastSquare )
      fResult->NormalizeErrors(); 

   
   
   // set also new parameter values and errors in FitConfig
   if (fConfig.UpdateAfterFit() && ret) DoUpdateFitConfig(); 

   return ret; 
}

bool Fitter::DoMinimization(const BaseFunc & objFunc, const ROOT::Math::IMultiGenFunction * chi2func) { 
   // perform the minimization initializing the minimizer starting from a given obj function 

   // keep also a copy of FCN function and set this in minimizer so they will be managed together
   // (remember that cloned copy will still depends on data and model function pointers) 
   fObjFunction = std::auto_ptr<ROOT::Math::IMultiGenFunction> ( objFunc.Clone() ); 
   if (!DoInitMinimizer()) return false;
   return DoMinimization(chi2func);
}


void Fitter::DoUpdateFitConfig() { 
   // update the fit configuration after a fit using the obtained result 
   if (fResult->IsEmpty() || !fResult->IsValid() ) return; 
   for (unsigned int i = 0; i < fConfig.NPar(); ++i) { 
      ParameterSettings & par = fConfig.ParSettings(i); 
      par.SetValue( fResult->Value(i) );
      if (fResult->Error(i) > 0) par.SetStepSize( fResult->Error(i) );
   } 
}

int Fitter::GetNCallsFromFCN() { 
   // retrieve ncalls from the fit method functions
   int ncalls = 0;
   if (!fUseGradient) { 
      const ROOT::Math::FitMethodFunction * fcn = dynamic_cast<const ROOT::Math::FitMethodFunction *>(fObjFunction.get());
      if (fcn) ncalls = fcn->NCalls();
   }
   else { 
      const ROOT::Math::FitMethodGradFunction * fcn = dynamic_cast<const ROOT::Math::FitMethodGradFunction*>(fObjFunction.get());
      if (fcn) ncalls = fcn->NCalls();
   }
   return ncalls;
}


bool Fitter::ApplyWeightCorrection(const ROOT::Math::IMultiGenFunction & loglw2, bool minimizeW2L) { 
   // apply correction for weight square 
   // Compute Hessian of the loglikelihood function using the sum of the weight squared 
   // This method assumes: 
   // - a fit has been done before and a covariance matrix exists
   // - the objective function is a likelihood function and Likelihood::UseSumOfWeightSquare()
   //    has been called before 

   if (fMinimizer.get() == 0) { 
      MATH_ERROR_MSG("Fitter::ApplyWeightCorrection","Must perform first a fit before applying the correction");
      return false;
   }

   unsigned int n = loglw2.NDim();
   // correct errors for weight squared 
   std::vector<double> cov(n*n);
   bool ret = fMinimizer->GetCovMatrix(&cov[0] );
   if (!ret) {
      MATH_ERROR_MSG("Fitter::ApplyWeightCorrection","Previous fit has no valid Covariance matrix");
      return false;
   }
   // need to re-init the minimizer and set w2
   fObjFunction = std::auto_ptr<ROOT::Math::IMultiGenFunction> ( loglw2.Clone() ); 
   // need to re-initialize the minimizer for the changes applied in the 
   // objective functions 
   if (!DoInitMinimizer()) return false;

   //std::cout << "Running Hesse ..." << std::endl;

   // run eventually before a minimization
   // ignore its error
   if (minimizeW2L) fMinimizer->Minimize();
   // run Hesse on the log-likelihood build using sum of weight squared 
   ret = fMinimizer->Hesse();
   if (!ret) { 
      MATH_ERROR_MSG("Fitter::ApplyWeightCorrection","Error running Hesse on weight2 likelihood - cannot compute errors");
      return false;
   }

   if (fMinimizer->CovMatrixStatus() != 3) { 
      MATH_WARN_MSG("Fitter::ApplyWeightCorrection","Covariance matrix for weighted likelihood is not accurate, the errors may be not reliable");
      if (fMinimizer->CovMatrixStatus() == 2)
         MATH_WARN_MSG("Fitter::ApplyWeightCorrection","Covariance matrix for weighted likelihood was forced to be defined positive");
      if (fMinimizer->CovMatrixStatus() <= 0)
         // probably should have failed before 
         MATH_ERROR_MSG("Fitter::ApplyWeightCorrection","Covariance matrix for weighted likelihood is not valid !");
   }
   
   // std::vector<double> c(n*n);
   //       ret = fMinimizer->GetCovMatrix(&c2[0] );
   //       if (!ret) std::cout << "Error reading cov matrix " << fMinimizer->Status() << std::endl;
   //       TMatrixDSym cmat2(n,&c2[0]);
   //       std::cout << "Cov matrix of w2 " << std::endl;
   //       cmat2.Print();
   //       cmat2.Invert();
   //       std::cout << "Hessian of w2 " << std::endl;
   //       cmat2.Print();

   // get Hessian matrix from weight-square likelihood
   std::vector<double> hes(n*n);
   ret = fMinimizer->GetHessianMatrix(&hes[0] );
   if (!ret) {
      MATH_ERROR_MSG("Fitter::ApplyWeightCorrection","Error retrieving Hesse on weight2 likelihood - cannot compute errors");
      return false;
   }

   // for debug
   // std::cout << "Hessian W2 matrix " << std::endl;
   // for (unsigned int i = 0; i < n; ++i) {
   //    for (unsigned int j = 0; j < n; ++j) {
   //       std::cout << std::setw(12) << hes[i*n + j] << " , ";
   //    }
   //    std::cout << std::endl;
   // }

   // perform product of matvrix cov * hes * cov
   // since we do not want to add matrix dependence do product by hand 
   // first do  hes * cov
   std::vector<double> tmp(n*n);
   for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
         for (unsigned int k = 0; k < n; ++k) 
            tmp[i*n+j] += hes[i*n + k] * cov[k*n + j];
      }
   }
   // do multiplication now cov * tmp save result 
   std::vector<double> newCov(n*n);
   for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
         for (unsigned int k = 0; k < n; ++k) 
            newCov[i*n+j] += cov[i*n + k] * tmp[k*n + j];
      }
   }
   // update fit result with new corrected covariance matrix
   unsigned int k = 0;
   for (unsigned int i = 0; i < n; ++i) { 
      fResult->fErrors[i] = std::sqrt( newCov[i*(n+1)] );
      for (unsigned int j = 0; j <= i; ++j) 
         fResult->fCovMatrix[k++] = newCov[i *n + j];
   }
   //fResult->PrintCovMatrix(std::cout);

   return true; 
}


   } // end namespace Fit

} // end namespace ROOT

