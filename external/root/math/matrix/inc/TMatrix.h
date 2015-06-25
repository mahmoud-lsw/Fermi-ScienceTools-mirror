// @(#)root/matrix:$Id: TMatrix.h,v 1.1.1.1 2009/04/05 20:44:15 elwinter Exp $
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrix
#define ROOT_TMatrix

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrix                                                              //
//                                                                      //
//  Instantation of TMatrixT<Float_t>                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMatrixF
#include "TMatrixF.h"
#endif
typedef TMatrixT<Float_t> TMatrix;

#endif
