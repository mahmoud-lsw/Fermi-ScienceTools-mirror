/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooMPSentinel.h,v 1.1.1.1 2009/04/05 20:44:27 elwinter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MP_SENTINEL
#define ROO_MP_SENTINEL

#include "Rtypes.h"
#include "RooArgSet.h"
class RooRealMPFE ;

class RooMPSentinel {
public:

  RooMPSentinel() ;
  virtual ~RooMPSentinel() ;
 
protected:

  friend class RooRealMPFE ;
  void add(RooRealMPFE& mpfe) ;
  void remove(RooRealMPFE& mpfe) ;

  RooMPSentinel(const RooMPSentinel&) {
    // Default constructor
  }
  RooArgSet _mpfeSet ;
  
  ClassDef(RooMPSentinel,1) // Singleton class that terminate MP server processes when parent exists
};

#endif
