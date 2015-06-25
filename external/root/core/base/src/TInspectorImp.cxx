// @(#)root/base:$Id: TInspectorImp.cxx,v 1.1.1.1 2009/04/05 20:44:31 elwinter Exp $
// Author: Fons Rademakers   07/05/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TInspector Imp                                                             //
//                                                                            //
// ABC describing GUI independent object inspector (abstration mainly needed  //
// for Win32. On X11 systems it currently uses a standard TCanvas).           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TInspectorImp.h"

ClassImp(TInspectorImp)
