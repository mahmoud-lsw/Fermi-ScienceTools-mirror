// @(#)root/matrix:$Id: TMatrixDBasefwd.h,v 1.1.1.2 2010/04/07 18:31:28 elwinter Exp $
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixDBasefwd
#define ROOT_TMatrixDBasefwd

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixDBase                                                         //
//                                                                      //
//  Forward declaration of TMatrixTBase<Double_t>                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

template<class Element> class TMatrixTBase;
typedef TMatrixTBase<Double_t> TMatrixDBase;

#endif
