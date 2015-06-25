// @(#)root/graf2d:$Id: QuartzUtils.h,v 1.1.1.1 2013/01/23 16:01:47 areustle Exp $
// Author: Timur Pocheptsov, 11/06/2012

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_QuartzUtils
#define ROOT_QuartzUtils

#include <ApplicationServices/ApplicationServices.h>

namespace ROOT {
namespace Quartz {

//Scope guard class for CGContextRef.
class CGStateGuard {
public:
   CGStateGuard(CGContextRef ctx);
   ~CGStateGuard();
   
private:
   CGContextRef fCtx;
   
   CGStateGuard(const CGStateGuard &rhs);
   CGStateGuard &operator = (const CGStateGuard &rhs);
};

}
}

#endif
