/* @(#)root/gdml:$Id: LinkDef.h,v 1.1.1.2 2013/01/23 16:01:13 areustle Exp $ */
/* Authors: Ben Lloyd 09/11/06 ben.lloyd@cern.ch */

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//#pragma link C++ class GDMLEngine;
#pragma link C++ class TGDMLParse;
#pragma link C++ class TGDMLRefl;
#pragma link C++ class TGDMLWrite;

#endif
