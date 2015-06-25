//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// version: $Id: CValuePrinter.h,v 1.1.1.1 2013/01/23 16:01:27 areustle Exp $
// author:  Vassil Vassilev <vasil.georgiev.vasilev@cern.ch>
//------------------------------------------------------------------------------

#ifndef CLING_VALUEPRINTERC_H
#define CLING_VALUEPRINTERC_H
#ifdef __cplusplus
extern "C"
#endif
void cling_PrintValue(void* /*clang::Expr**/ E,
                      void* /*clang::ASTContext**/ C,
                      const void* value);

#endif // CLING_VALUEPRINTERC_H
