// @(#)root/graf2d:$Id: QuartzMarker.h,v 1.1.1.1 2013/01/23 16:01:47 areustle Exp $
// Author: Timur Pocheptsov, 14/8/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_QuartzMarker
#define ROOT_QuartzMarker

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// QuartzMarkers                                                        //
//                                                                      //
// Aux. functions to draw poly-markers.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include <Cocoa/Cocoa.h>

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#ifndef ROOT_TPoint
#include "TPoint.h"
#endif

namespace ROOT {
namespace Quartz {

void DrawPolyMarker(CGContextRef ctx, const std::vector<TPoint> &marker, 
                    Size_t markerSize, Style_t markerStyle);
void DrawPolyMarker(CGContextRef ctx, unsigned nPoints, const TPoint *marker, 
                    Size_t markerSize, Style_t markerStyle);

}
}

#endif
