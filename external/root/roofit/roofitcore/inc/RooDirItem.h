/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDirItem.h,v 1.1.1.1 2009/04/05 20:44:27 elwinter Exp $
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
#ifndef ROO_DIR_ITEM
#define ROO_DIR_ITEM

#include "Rtypes.h"
class TDirectory ;

class RooDirItem {
public:
  RooDirItem() ;
  RooDirItem(const RooDirItem& other) ;
  virtual ~RooDirItem() ;

protected:

  void appendToDir(TObject* obj, Bool_t forceMemoryResident=kFALSE) ;
  void removeFromDir(TObject* obj) ;

  TDirectory* _dir ;     //! Associated directory
  ClassDef(RooDirItem,1) // Base class for RooFit objects that are listed TDirectories
};

#endif
