// @(#)root/treeplayer:$Id: TTreeProxyGenerator.h,v 1.1.1.2 2013/01/23 16:01:06 areustle Exp $
// Author: Philippe Canal 01/06/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeProxyGenerator
#define ROOT_TTreeProxyGenerator

#ifndef ROOT_Tlist
#include "TList.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

class TBranch;
class TBranchElement;
class TLeaf;
class TTree;
class TStreamerElement;

namespace ROOT {
   class TBranchProxy;
   class TFriendProxyDescriptor;
   class TBranchProxyDescriptor;
   class TBranchProxyClassDescriptor;

   class TTreeProxyGenerator
   {
   public:
      enum EContainer { kNone, kClones, kSTL };
      enum EOption { kNoOption, kNoHist };
      UInt_t   fMaxDatamemberType;
      TString  fScript;
      TString  fCutScript;
      TString  fPrefix;
      TString  fHeaderFileName;
      TString  fOptionStr;
      UInt_t   fOptions;
      UInt_t   fMaxUnrolling;
      TTree   *fTree;
      TList    fListOfHeaders;
      TList    fListOfClasses;
      TList    fListOfFriends;
      TList    fListOfPragmas;
      TList    fListOfTopProxies;
      TList   *fCurrentListOfTopProxies; //!
      TList    fListOfForwards;
      TTreeProxyGenerator(TTree* tree, const char *script, const char *fileprefix, 
                          const char *option, UInt_t maxUnrolling);
      TTreeProxyGenerator(TTree* tree, const char *script, const char *cutscript, 
                          const char *fileprefix, const char *option, UInt_t maxUnrolling);

      TBranchProxyClassDescriptor* AddClass(TBranchProxyClassDescriptor *desc);
      void AddDescriptor(TBranchProxyDescriptor *desc);
      void AddForward(TClass *cl);
      void AddForward(const char *classname);
      void AddFriend(TFriendProxyDescriptor *desc);
      void AddHeader(TClass *cl);
      void AddHeader(const char *classname);
      void AddMissingClassAsEnum(const char *clname, Bool_t isscope);
      void AddPragma(const char *pragma_text);
      void CheckForMissingClass(const char *clname);

      Bool_t NeedToEmulate(TClass *cl, UInt_t level);

      void   ParseOptions();

      UInt_t AnalyzeBranches(UInt_t level, TBranchProxyClassDescriptor *topdesc, TBranchElement *branch, TVirtualStreamerInfo *info = 0);
      UInt_t AnalyzeBranches(UInt_t level, TBranchProxyClassDescriptor *topdesc, TIter &branches, TVirtualStreamerInfo *info);
      UInt_t AnalyzeOldBranch(TBranch *branch, UInt_t level, TBranchProxyClassDescriptor *desc);
      UInt_t AnalyzeOldLeaf(TLeaf *leaf, UInt_t level, TBranchProxyClassDescriptor *topdesc);
      void   AnalyzeElement(TBranch *branch, TStreamerElement *element, UInt_t level, TBranchProxyClassDescriptor *desc, const char* path);
      void   AnalyzeTree(TTree *tree);
      void   WriteProxy();

      const char *GetFileName() { return fHeaderFileName; }
   };

}

using ROOT::TTreeProxyGenerator;

#endif