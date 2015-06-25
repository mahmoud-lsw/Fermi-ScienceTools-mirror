/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooTable.h,v 1.1.1.1 2009/04/05 20:44:26 elwinter Exp $
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
#ifndef ROO_TABLE
#define ROO_TABLE

#include "Riosfwd.h"
#include <assert.h>
#include "TNamed.h"
#include "RooAbsCategory.h"
#include "RooPrintable.h"

class RooTable : public TNamed, public RooPrintable {
public:

  // Constructors, cloning and assignment
  RooTable() {} ;
  virtual ~RooTable() ;
  RooTable(const char *name, const char *title);
  RooTable(const RooTable& other) ;

  virtual void fill(RooAbsCategory& cat, Double_t weight=1.0) = 0 ;

  virtual Bool_t isIdentical(const RooTable& other) = 0 ;

protected:

  ClassDef(RooTable,1) // Abstract interface for tables
};

#endif