// -*- mode: c++ -*-
/**
 * @file pyExposure.in
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /glast/ScienceTools/glast/pyExposure/src/pyExposure.i,v 1.1.1.2 2011/03/20 19:24:45 elwinter Exp $
 */
%module pyExposure
%{
#ifdef TRAP_FPE
#include <fenv.h>
#endif
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "evtbin/Gti.h"
#include "dataSubselector/Cuts.h"
#include "dataSubselector/GtiCut.h"
#include "dataSubselector/RangeCut.h"
#include "dataSubselector/SkyConeCut.h"
#include "pyExposure/Exposure.h"
%}
%include stl.i
%exception {
   try {
      $action
   } catch (std::exception & eObj) {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(eObj.what()));
      return NULL;
   }
}
%template(DoubleVector) std::vector<double>;
%template(DoubleVectorVector) std::vector< std::vector<double> >;
%template(StringVector) std::vector<std::string>;
%include evtbin/Gti.h
%rename (evtbin_Gti) Gti;
%include dataSubselector/CutBase.h
%include dataSubselector/Gti.h
%include dataSubselector/GtiCut.h
%include dataSubselector/RangeCut.h
%include dataSubselector/SkyConeCut.h
%include dataSubselector/Cuts.h
%template(GtiCutVector) std::vector<dataSubselector::GtiCut *>;
%include pyExposure/Exposure.h
#ifdef TRAP_FPE
%extend pyExposure::Exposure {
   static void enableFPE() {
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
   }
}
#endif
%extend dataSubselector::Cuts {
   void writeCuts() const {
      self->writeCuts(std::cout);
   }
   std::vector<dataSubselector::GtiCut *> gtiCuts() {
      std::vector<const dataSubselector::GtiCut *> my_gtiCuts;
      self->getGtiCuts(my_gtiCuts);
      std::vector<dataSubselector::GtiCut *> output;
      for (size_t i = 0; i < my_gtiCuts.size(); i++) {
         output.push_back(my_gtiCuts.at(i)->clone());
      }
      return output;
   }
   dataSubselector::CutBase & getCut(size_t i) {
      return const_cast<dataSubselector::CutBase &>(self->operator[](i));
   }
   static dataSubselector::SkyConeCut * 
   castAsSkyConeCut(dataSubselector::CutBase * cut) {
      dataSubselector::SkyConeCut * my_cut 
         = dynamic_cast<dataSubselector::SkyConeCut *>(cut);
      if (cut == 0) {
         throw std::runtime_error("Cannot cast to SkyConeCut.");
      }
      return my_cut;
   }
   static dataSubselector::RangeCut * 
   castAsRangeCut(dataSubselector::CutBase * cut) {
      dataSubselector::RangeCut * my_cut 
         = dynamic_cast<dataSubselector::RangeCut *>(cut);
      if (cut == 0) {
         throw std::runtime_error("Cannot cast to RangeCut.");
      }
      return my_cut;
   }
}

