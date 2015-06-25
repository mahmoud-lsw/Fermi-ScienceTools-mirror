//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// version: $Id: ASTUtils.h,v 1.1.1.1 2013/01/23 16:01:27 areustle Exp $
// author:  Vassil Vassilev <vasil.georgiev.vasilev@cern.ch>
//------------------------------------------------------------------------------

#ifndef CLING_AST_UTILS_H
#define CLING_AST_UTILS_H

#include "clang/AST/Type.h"

namespace clang {
  class Expr;
  class DeclContext;
  class NamedDecl;
  class NamespaceDecl;
  class Sema;
}

namespace cling {
  class Synthesize {
  public:
    static clang::Expr* CStyleCastPtrExpr(clang::Sema* S, 
                                          clang::QualType Ty, uint64_t Ptr);
  };

  class Lookup {
  public:
    static clang::NamespaceDecl* Namespace(clang::Sema* S, 
                                           const char* Name,
                                           clang::DeclContext* Within = 0);
    static clang::NamedDecl* Named(clang::Sema* S, 
                                   const char* Name,
                                   clang::DeclContext* Within = 0);

  };
} // end namespace cling
#endif // CLING_AST_UTILS_H
