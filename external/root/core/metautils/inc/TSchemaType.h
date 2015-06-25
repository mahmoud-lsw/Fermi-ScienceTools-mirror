// @(#)root/core:$Id: TSchemaType.h,v 1.1.1.1 2013/02/07 21:31:30 areustle Exp $

#ifndef R__TSCHEMATYPE_H
#define R__TSCHEMATYPE_H

#if !defined(__CINT__)
// Avoid clutering the dictionary (in particular with the STL declaration)

#include <string>

namespace ROOT
{
   struct TSchemaType {
      TSchemaType() {}
      TSchemaType(const char *type, const char *dim) : fType(type),fDimensions(dim) {}
      TSchemaType(std::string &type, std::string &dim) : fType(type),fDimensions(dim) {}
      std::string fType;
      std::string fDimensions;
   };
}
#endif // __CINT__
#endif // R__TSCHEMATYPE_H

