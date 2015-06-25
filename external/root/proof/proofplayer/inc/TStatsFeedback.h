// @(#)root/proofplayer:$Id: TStatsFeedback.h,v 1.1.1.1 2013/01/23 16:01:31 areustle Exp $
// Author: G. Ganis May 2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStatsFeedback
#define ROOT_TStatsFeedback


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStatsFeedback                                                       //
//                                                                      //
// Utility class to display PROOF stats feedback histos during queries. //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TQObject
#include "TQObject.h"
#endif
#ifndef ROOT_TObject
#include "TObject.h"
#endif

class TProof;
class TProof;
class THashList;
class TSeqCollection;

class TStatsFeedback : public TObject, public TQObject {
protected:
   TString        fName;   // Identfier for this object
   TProof        *fProof;  //handle to PROOF session
public:
   TStatsFeedback(TProof *proof = 0);
   ~TStatsFeedback();

   void        Feedback(TList *objs);
   const char *GetName() const { return fName.Data(); }
   ULong_t     Hash() const { return fName.Hash(); }

   ClassDef(TStatsFeedback,0)  // Present PROOF query feedback
};

#endif