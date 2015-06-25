//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// version: $Id: InvocationOptions.h,v 1.1.1.1 2013/01/23 16:01:27 areustle Exp $
// author:  Axel Naumann <axel@cern.ch>
//------------------------------------------------------------------------------

#ifndef CLING_INVOCATIONOPTIONS_H
#define CLING_INVOCATIONOPTIONS_H

#include "llvm/Support/Path.h"

#include <vector>

namespace cling {
  class InvocationOptions {
  public:
    InvocationOptions():
      NoLogo(false), ShowVersion(false), Verbose(false), Help(false) {}
    bool NoLogo;
    bool ShowVersion;
    bool Verbose;
    bool Help;

    std::vector<std::string> LibsToLoad;
    std::vector<llvm::sys::Path> LibSearchPath;

    static InvocationOptions CreateFromArgs(int argc, const char* const argv[],
                                            std::vector<unsigned>& leftoverArgs
                                            /* , Diagnostic &Diags */);

    void PrintHelp();
  };
}

#endif // INVOCATIONOPTIONS
