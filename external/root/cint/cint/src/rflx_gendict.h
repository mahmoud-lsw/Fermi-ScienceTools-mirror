/* -*- C++ -*- */
/*************************************************************************
 * Copyright(c) 1995~2005  Masaharu Goto (root-cint@cern.ch)
 *
 * For the licensing terms see the file COPYING
 *
 ************************************************************************/
//$Id: rflx_gendict.h,v 1.1.1.2 2013/01/23 16:01:19 areustle Exp $

#ifndef RFLX_GENDICT_H
#define RFLX_GENDICT_H 1

#ifdef __cplusplus
extern "C" {
#endif 

  void rflx_gendict(const char * linkfilename,
		    const char * sourcefile);

#ifdef __cplusplus
}
#endif 

#endif // RFLX_GENDICT_H
