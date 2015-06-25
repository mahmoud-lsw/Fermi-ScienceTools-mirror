// @(#)root/base:$Id: TVirtualMonitoring.cxx,v 1.1.1.1 2009/04/05 20:44:31 elwinter Exp $
// Author: Andreas-Joachim Peters  15/05/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualMonitoring                                                   //
//                                                                      //
// Provides the interface for externel Monitoring                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TVirtualMonitoring.h"


ClassImp(TVirtualMonitoringWriter)
ClassImp(TVirtualMonitoringReader)


TVirtualMonitoringWriter *gMonitoringWriter = 0;
TVirtualMonitoringReader *gMonitoringReader = 0;

