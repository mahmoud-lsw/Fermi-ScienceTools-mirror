// @(#)root/metautils:$Id: TMetaUtils.h,v 1.1.1.1 2013/02/07 21:31:30 areustle Exp $
// Author: Axel Naumann, Nov 2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMetaUtils
#define ROOT_TMetaUtils

// All clang entities must stay opaque types
namespace clang {
   class QualType;
}
namespace cling {
   class Interpreter;
}

namespace ROOT {
   namespace TMetaUtils {
      clang::QualType LookupTypeDecl(cling::Interpreter& interp,
                                     const char* tyname);
   }
}
#endif // ROOT_TMetaUtils
