// @(#)root/net:$Id: TGSFile.cxx,v 1.1.1.1 2013/01/23 16:01:06 areustle Exp $
// Author: Marcelo Sousa   23/08/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGSFile                                                              //
//                                                                      //
// A TGSFile is a normal TWebFile but it reads data from the            //
// Google Storage server. As a derived TWebFile class TGSFile it is     //
// a read only file. The HTTP requests are generated by THTTPMessage    //
// objects with the auth_prefix set as GOOG1. The user id and secret    //
// pass required to sign the requests are passed through the            //
// environment variables GT_ACCESS_ID and GT_ACCESS_KEY.                //
// For more information check:                                          //
//   http://code.google.com/apis/storage/docs/getting-started.html      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGSFile.h"
#include "THTTPMessage.h"
#include "TSocket.h"
#include "TROOT.h"

#include <errno.h>
#include <stdlib.h>

ClassImp(TGSFile)

//______________________________________________________________________________
TGSFile::TGSFile(const char *path, Option_t *) : TWebFile(path, "IO")
{
   // For TGSFile to properly work you need to set up
   // environment variables GT_ACCESS_ID and GT_ACCESS_KEY.
   // The format of the path is: server/bucket/file, e.g.
   // f = new TGSFile("gs://commondatastorage.googleapis.com/roots3/hsimple.root")

   TString tpath = path;

   Int_t begPath = 5, slash = 0, i = 0;

   if (tpath.BeginsWith("gs://") == kFALSE) {
      Error("TGSFile", "invalid path %s", path);
      goto zombie;
   }

   while (i < 2 && begPath < tpath.Length()) {
      slash = tpath.Index('/', begPath);

      if (slash == kNPOS) {
         Error("TGSFile","invalid path %s", path);
         goto zombie;
      }

      switch(i){
         case 0:
            fServer = TUrl(TString(tpath(begPath,slash)));
            break;
         case 1:
            fBucket = tpath(begPath,slash-begPath);
            fRealName = "/" + tpath(slash+1, tpath.Length()-(slash+1));
      }
      i++;
      begPath = slash+1;
   }

   fAuthPrefix = "GOOG1";
   fAccessId   = gSystem->Getenv("GT_ACCESS_ID");
   fAccessKey  = gSystem->Getenv("GT_ACCESS_KEY");
   if (fAccessId == "" || fAccessKey == "") {
      if (fAccessId == "")  Error("TGSFile", "shell variable GT_ACCESS_ID not set");
      if (fAccessKey == "") Error("TGSFile", "shell variable GT_ACCESS_KEY not set");
      goto zombie;
   }

   Init(kFALSE);
   return;
   
zombie:
   MakeZombie();
   gDirectory = gROOT;
}

//______________________________________________________________________________
Int_t TGSFile::GetHead()
{
   // Clone of TWebFile::GetHead except it uses THTTPMessage to generate
   // the HTTP request.

   THTTPMessage s3head = THTTPMessage(kHEAD, fRealName, GetBucket(),
                                      GetUrl().GetHost(), GetAuthPrefix(),
                                      GetAccessId(), GetAccessKey());

   TString msg = s3head.GetRequest();

   TUrl connurl;

   fUrl = fServer;

   if (fProxy.IsValid())
      connurl = fProxy;
   else
      connurl = fUrl;

   TSocket *s = 0;
   for (Int_t i = 0; i < 5; i++) {
      s = new TSocket(connurl.GetHost(), connurl.GetPort());
      if (!s->IsValid()) {
         delete s;
         if (gSystem->GetErrno() == EADDRINUSE || gSystem->GetErrno() == EISCONN) {
            s = 0;
            gSystem->Sleep(i*10);
         } else {
            Error("GetHead", "cannot connect to host %s (errno=%d)", fUrl.GetHost(),
                  gSystem->GetErrno());
            return -1;
         }
      } else
         break;
   }
   if (!s)
      return -1;

   if (s->SendRaw(msg.Data(), msg.Length()) == -1) {
      Error("GetHead", "error sending command to host %s", fUrl.GetHost());
      delete s;
      return -1;
   }

   char line[8192];
   Int_t n, ret = 0, redirect = 0;

   while ((n = GetLine(s, line, sizeof(line))) >= 0) {
      if (n == 0) {
         if (gDebug > 0)
            Info("GetHead", "got all headers");
         delete s;
         if (fBasicUrlOrg != "" && !redirect) {
            // set back to original url in case of temp redirect
            SetMsgReadBuffer10();
            fMsgGetHead = "";
         }
         if (ret < 0)
            return ret;
         if (redirect)
            return GetHead();
         return 0;
      }

      if (gDebug > 0)
         Info("GetHead", "header: %s", line);

      TString res = line;
      if (res.BeginsWith("HTTP/1.")) {
         if (res.BeginsWith("HTTP/1.1")) {
            if (!fHTTP11) {
               fMsgGetHead = "";
               fMsgReadBuffer10 = "";
            }
            fHTTP11 = kTRUE;
         }
         TString scode = res(9, 3);
         Int_t code = scode.Atoi();
         if (code >= 500) {
            if (code == 500)
               fHasModRoot = kTRUE;
            else {
               ret = -1;
               TString mess = res(13, 1000);
               Error("GetHead", "%s: %s (%d)", fBasicUrl.Data(), mess.Data(), code);
            }
         } else if (code >= 400) {
            if (code == 400)
               ret = -3;   // command not supported
            else if (code == 404)
               ret = -2;   // file does not exist
            else {
               ret = -1;
               TString mess = res(13, 1000);
               Error("GetHead", "%s: %s (%d)", fBasicUrl.Data(), mess.Data(), code);
            }
         } else if (code >= 300) {
            if (code == 301 || code == 303)
               redirect = 1;   // permanent redirect
            else if (code == 302 || code == 307)
               redirect = 2;   // temp redirect
            else {
               ret = -1;
               TString mess = res(13, 1000);
               Error("GetHead", "%s: %s (%d)", fBasicUrl.Data(), mess.Data(), code);
            }
         } else if (code > 200) {
            ret = -1;
            TString mess = res(13, 1000);
            Error("GetHead", "%s: %s (%d)", fBasicUrl.Data(), mess.Data(), code);
         }
      } else if (res.BeginsWith("Content-Length:")) {
         TString slen = res(16, 1000);
         fSize = slen.Atoll();
      } else if (res.BeginsWith("Location:") && redirect) {
         TString redir = res(10, 1000);
         if (redirect == 2)   // temp redirect
            SetMsgReadBuffer10(redir, kTRUE);
         else               // permanent redirect
            SetMsgReadBuffer10(redir, kFALSE);
         fMsgGetHead = "";
      }
   }

   delete s;
   return ret;
}

#if 0
// Currently not supported on the Google Storage cloud side because
// the chunked responses are not "standard"
//______________________________________________________________________________
Bool_t  TGSFile::ReadBuffers10(char *buf, Long64_t *pos, Int_t *len, Int_t nbuf)
{
   // Read specified byte range from Google Storage.
   // This routine connects to the Google Storage server, sends the
   // request created by THTTPMessage and returns the buffer.
   // Returns kTRUE in case of error.
	
   THTTPMessage gsget = THTTPMessage(kGET, fRealName, GetBucket(),
                                     GetUrl().GetHost(), GetAuthPrefix(),
                                     GetAccessId(), GetAccessKey(),
                                     fOffset, pos, len, nbuf);
   TString msg = gsget.GetRequest();
   Int_t size = gsget.GetLength();
   Int_t curbuf = gsget.GetCurrentBuffer();
   
   //printf("Num of buffers %d, current buffer %d size %d\n",nbuf,curbuf,size);
   
   Int_t n = GetFromWeb10(buf, size, msg);
   if (n == -1)
      return kTRUE;
   // The -2 error condition typically only happens when
   // GetHead() failed because not implemented, in the first call to
   // ReadBuffer() in Init(), it is not checked in ReadBuffers10().
   if (n == -2) {
      Error("ReadBuffer10", "%s does not exist", fBasicUrl.Data());
      MakeZombie();
      gDirectory = gROOT;
      return kTRUE;
   }

   fOffset += size;

   if (nbuf == curbuf) {
      return kFALSE;         
   } else {
      return ReadBuffers10(&buf[size], &pos[curbuf], &len[curbuf], nbuf-curbuf);
   }   
}
#endif

//______________________________________________________________________________
Bool_t TGSFile::ReadBuffer(char *buf, Int_t len)
{
   return ReadBuffer10(buf,len);
}

//______________________________________________________________________________
Bool_t TGSFile::ReadBuffer10(char *buf, Int_t len)
{
   // Read specified byte range from Google Storage.
   // This routine connects to the Google Storage server, sends the
   // request created by THTTPMessage and returns the buffer.
   // Returns kTRUE in case of error.

   const Int_t nbuf = 1;
   Long64_t pos[nbuf];
   pos[nbuf-1] = fOffset; 

   Int_t lens[nbuf];
   lens[nbuf-1] = len;

   THTTPMessage gsget = THTTPMessage(kGET, fRealName, GetBucket(),
                                     GetUrl().GetHost(), GetAuthPrefix(),
                                     GetAccessId(), GetAccessKey(),
                                     0, pos, lens, nbuf);
   TString msg = gsget.GetRequest();

   Int_t n = GetFromWeb10(buf, len, msg);
   if (n == -1)
      return kTRUE;
   // The -2 error condition typically only happens when
   // GetHead() failed because not implemented, in the first call to
   // ReadBuffer() in Init(), it is not checked in ReadBuffers10().
   if (n == -2) {
      Error("ReadBuffer10", "%s does not exist", fBasicUrl.Data());
      MakeZombie();
      gDirectory = gROOT;
      return kTRUE;
   }

   fOffset += len;

   return kFALSE;
}

