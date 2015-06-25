/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooFit.h,v 1.1.1.1 2009/04/05 20:44:27 elwinter Exp $
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
#ifndef ROO_FIT
#define ROO_FIT

#include "RConfig.h"

// Global include file to fix occasional compiler issues
// An error in the construction of the system and C++ header files on
// Solaris 8 / Workshop 6 Updates 1&2 leads to a conflict between the use
// of ::clock_t and clock_t when <string> is compiled under
// -D_XOPEN_SOURCE=500. The following code ensures that ::clock_t is
// always defined and thus allows <string> to compile.
// This is just a workaround and should be monitored as compiler and
// operating system versions evolve.
#if defined(__SUNPRO_CC) && defined(_XOPEN_SOURCE) && (_XOPEN_SOURCE - 0 == 500 )
#ifndef _CLOCK_T
#define _CLOCK_T
typedef long clock_t; /* relative time in a
specified resolution */
#endif /* ifndef _CLOCK_T */
#endif // SUN and XOPENSOURCE=500


#endif
