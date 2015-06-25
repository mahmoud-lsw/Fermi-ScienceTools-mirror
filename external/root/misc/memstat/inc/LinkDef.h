// @(#)root/memstat:$Name: ScienceTools-v10r0p5-fssc-20150518 $:$Id: LinkDef.h,v 1.1.1.2 2013/01/23 16:01:12 areustle Exp $
// Author: Anar Manafov (A.Manafov@gsi.de) 2008-04-28

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace memstat;

#pragma link C++ class TMemStat;
#pragma link C++ class memstat::TMemStatMng;


//#pragma link C++ function Memstat::dig2bytes(Long64_t);

#endif
