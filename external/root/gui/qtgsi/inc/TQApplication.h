// @(#)root/qtgsi:$Id: TQApplication.h,v 1.1.1.1 2009/04/05 20:43:53 elwinter Exp $
// Author: Denis Bertini, M. Al-Turany  01/11/2000

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TQApplication
#define ROOT_TQApplication

//////////////////////////////////////////////////////////////////////////////////
//
//  TQApplication
//
//  This class creates ROOT environement that will
//  interface with the Qt windowing system eventloop and eventhandlers.
//  This class will be instantiated once (singleton) in a main()
//  program.
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TApplication
#include "TApplication.h"
#endif

class TQApplication : public TApplication {

private:
   Bool_t fCustomized; // flag for customized canvas implimentation

public:
   TQApplication();
   TQApplication(const char *appClassName, int *argc, char **argv, void *options = 0, int numOptions = 0);
   virtual ~TQApplication();
   virtual void LoadGraphicsLibs();
   void SetCustomized();
   
   ClassDef(TQApplication,0) //creates ROOT environement with the Qt windowing system
};

#endif