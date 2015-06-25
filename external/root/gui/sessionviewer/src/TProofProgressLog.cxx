// @(#)root/sessionviewer:$Id: TProofProgressLog.cxx,v 1.1.1.3 2013/01/23 16:00:49 areustle Exp $
// Author: G Ganis, Jul 2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TError.h"
#include "TGFrame.h"
#include "TGTextView.h"
#include "TGScrollBar.h"
#include "TGLabel.h"
#include "TProof.h"
#include "TProofProgressDialog.h"
#include "TProofProgressLog.h"
#include "TProofLog.h"
#include "TGNumberEntry.h"
#include "TGListBox.h"
#include "TGMenu.h"
#include "TGButton.h"

const UInt_t kLogElemFilled = BIT(17); // If the log element has been retrieved at least once
const UInt_t kDefaultActive = BIT(18); // If the log element is active by default

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TProofProgressLog                                                     //
//                                                                       //
// Dialog used to display Proof session logs from the Proof progress     //
// dialog.                                                               //
// It uses TProofMgr::GetSessionLogs() mechanism internally              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(TProofProgressLog)

//____________________________________________________________________________
TProofProgressLog::TProofProgressLog(TProofProgressDialog *d, Int_t w, Int_t h) :
   TGTransientFrame(gClient->GetRoot(), gClient->GetRoot(), w, h)
{
   // Create a window frame for log messages.

   fDialog = d;
   if (fDialog) fSessionUrl = fDialog->fSessionUrl;
   fSessionIdx = 0;

   Init(w, h);
}

//____________________________________________________________________________
TProofProgressLog::TProofProgressLog(const char *url, Int_t idx, Int_t w, Int_t h) :
   TGTransientFrame(gClient->GetRoot(), gClient->GetRoot(), w, h)
{
   // Create a window frame for log messages.

   fDialog = 0;
   fSessionUrl = url;
   fSessionIdx = (idx > 0) ? -idx : idx;

   Init(w, h);
}

//____________________________________________________________________________
void TProofProgressLog::Init(Int_t w, Int_t h)
{
   // Init window frame for log messages.

   fProofLog = 0;
   fFullText = kFALSE;
   fTextType = kStd;
   // use hierarchical cleaning
   SetCleanup(kDeepCleanup);

   //The text window
   TGHorizontalFrame *htotal = new TGHorizontalFrame(this, w, h);
   TGVerticalFrame *vtextbox = new TGVerticalFrame(htotal, w, h);
   //fText = new TGTextView(this, w, h);
   fText = new TGTextView(vtextbox, w, h);
   vtextbox->AddFrame(fText, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3));

   //The frame for choosing workers
   fVworkers = new TGVerticalFrame(htotal);
   // URL choice
   TGLabel *laburl = new TGLabel(fVworkers, "Enter cluster URL:");
   fVworkers->AddFrame(laburl, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 2, 2, 2));
   fUrlText = new TGTextEntry(fVworkers);
   fUrlText->SetText(fSessionUrl.Data());
   fVworkers->AddFrame(fUrlText, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 0, 0, 0));
   //The lower row of number entries and buttons
   TGHorizontalFrame *hfurlbox = new TGHorizontalFrame(fVworkers, 20, 20);
   TGLabel *labsess = new TGLabel(hfurlbox, "Enter session:");
   hfurlbox->AddFrame(labsess, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 2, 2, 2));
   fSessNum = new TGNumberEntry(hfurlbox, 0, 5, -1, TGNumberFormat::kNESInteger);
   fSessNum->SetLimits(TGNumberFormat::kNELLimitMax, 0., 0.);
   fSessNum->SetIntNumber(0);
   fSessNum->GetNumberEntry()->SetToolTipText("Use 0 for the last known one,"
                                              " negative numbers for the previous ones, e.g. -1 for the last-but-one");
   hfurlbox->AddFrame(fSessNum, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 0, 0, 0));
   fUrlButton = new TGTextButton(hfurlbox, "Get logs info");
   fUrlButton->Connect("Clicked()", "TProofProgressLog", this, "Rebuild()");
   hfurlbox->AddFrame(fUrlButton, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 4, 0, 0, 0));
   fVworkers->AddFrame(hfurlbox, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));

   TGNumberEntry *nent = new TGNumberEntry(hfurlbox);
   fVworkers->AddFrame(nent, new TGLayoutHints(kLHintsTop | kLHintsLeft, 4, 0, 0, 0));

   //The list of workers
   fLogList = 0;
   BuildLogList(kTRUE);
   fLogList->Resize(102,52);
   fLogList->SetMultipleSelections(kTRUE); 

   //The SelectAll/ClearAll buttons
   TGHorizontalFrame *hfselbox = new TGHorizontalFrame(fVworkers, 20, 20);
   TGLabel *label1 = new TGLabel(hfselbox,"Choose workers:");
   hfselbox->AddFrame(label1, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 0, 0, 0));
   TGTextButton *selall = new TGTextButton(hfselbox,   "     &All      ");
   selall->Connect("Clicked()", "TProofProgressLog", this, "Select(=0)");
   hfselbox->AddFrame(selall, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 10, 0, 0, 0));
   TGTextButton *clearall = new TGTextButton(hfselbox, "     &Clear    ");
   clearall->Connect("Clicked()", "TProofProgressLog", this, "Select(=1)");
   hfselbox->AddFrame(clearall, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 10, 0, 0, 0));

   //select the defaut actives to start with
   Select(0, kFALSE);

   //Display button
   fLogNew = new TGTextButton(fVworkers, "&Display");
   fLogNew->Connect("Clicked()", "TProofProgressLog", this, "DoLog(=kFALSE)");
   //fLogNew->Resize(102, 20);
   // fLogNew->SetMargins(1, 1, 0, 1);
   fLogNew->SetTextColor(0xffffff, kFALSE);
   fLogNew->SetBackgroundColor(0x000044);
   fVworkers->AddFrame(hfselbox, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5, 2, 2, 2));
   fVworkers->AddFrame(fLogList, new TGLayoutHints(kLHintsExpandX | kLHintsTop | kLHintsExpandY, 2, 2, 5, 2));
   fVworkers->AddFrame(fLogNew, new TGLayoutHints(kLHintsExpandX | kLHintsTop , 2, 2, 1, 5));

   htotal->AddFrame(fVworkers, new TGLayoutHints(kLHintsCenterY | kLHintsLeft | kLHintsExpandY, 2, 2, 2, 2));

   //The lower row of number entries and buttons
   TGHorizontalFrame *hflogbox = new TGHorizontalFrame(vtextbox, 550, 20);
   fClose = new TGTextButton(hflogbox, "  &Close  ");
   fClose->Connect("Clicked()", "TProofProgressLog", this, "CloseWindow()");
   hflogbox->AddFrame(fClose, new TGLayoutHints(kLHintsCenterY |
                                                kLHintsRight, 10, 2, 2, 2));

   //Saving to a file controls
   fSave = new TGTextButton(hflogbox, "&Save");
   fSave->Connect("Clicked()", "TProofProgressLog", this, "SaveToFile()");
   hflogbox->AddFrame(fSave, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 4, 0, 0, 0));
   fFileName = new TGTextEntry(hflogbox);
   fFileName->SetText("<session-tag>.log");
   hflogbox->AddFrame(fFileName, new TGLayoutHints(kLHintsCenterY | kLHintsRight));
   TGLabel *label10 = new TGLabel(hflogbox, "Save to a file:");
   hflogbox->AddFrame(label10, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 50, 2, 2, 2));

   //Choose the number of lines to display
   TGVerticalFrame *vlines = new TGVerticalFrame(hflogbox);
   TGHorizontalFrame *vlines_buttons = new TGHorizontalFrame(vlines);
   TGLabel *label2 = new TGLabel(vlines_buttons, "Lines:");
   vlines_buttons->AddFrame(label2, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   fAllLines = new TGCheckButton(vlines_buttons, "all");
   fAllLines->SetToolTipText("Retrieve all lines (service messages excluded)");
   fAllLines->SetState(kButtonUp);
   fAllLines->Connect("Clicked()", "TProofProgressLog", this, "NoLineEntry()");
   vlines_buttons->AddFrame(fAllLines, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   fRawLines = new TGCheckButton(vlines_buttons, "svcmsg");
   fRawLines->SetToolTipText("Retrieve all type of lines, service messages included");
   fRawLines->SetState(kButtonUp);
   vlines_buttons->AddFrame(fRawLines, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   TGLabel *label11 = new TGLabel(vlines_buttons, "From");
   vlines_buttons->AddFrame(label11, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   fLinesFrom = new TGNumberEntry(vlines_buttons, 0, 5, -1, TGNumberFormat::kNESInteger);
   // coverity[negative_returns]: no problem with -100, the format is kNESInteger
   fLinesFrom->SetIntNumber(-100);
   fLinesFrom->GetNumberEntry()->SetToolTipText("Negative values indicate \"tail\" action");


   vlines_buttons->AddFrame(fLinesFrom, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   TGLabel *label3 = new TGLabel(vlines_buttons, "to");
   vlines_buttons->AddFrame(label3, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));
   fLinesTo = new TGNumberEntry(vlines_buttons, 0, 5, -1, TGNumberFormat::kNESInteger);
   vlines_buttons->AddFrame(fLinesTo, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));
   vlines->AddFrame(vlines_buttons, new TGLayoutHints(kLHintsCenterY));
   hflogbox->AddFrame(vlines, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

   //Grep controls
   TGLabel *label4 = new TGLabel(hflogbox, "Grep for:");
   hflogbox->AddFrame(label4, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 2, 2, 2));
   fGrepText = new TGTextEntry(hflogbox);
   hflogbox->AddFrame(fGrepText, new TGLayoutHints(kLHintsCenterY | kLHintsLeft));

   fGrepButton = new TGTextButton(hflogbox, "Grep");
   fGrepButton->Connect("Clicked()", "TProofProgressLog", this, "DoLog(=kTRUE)");
   hflogbox->AddFrame(fGrepButton, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 4, 0, 0, 0));

   vtextbox->AddFrame(hflogbox, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
   htotal->AddFrame(vtextbox, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY | kLHintsRight, 3, 3, 3, 3));
   AddFrame(htotal, new TGLayoutHints(kLHintsExpandX |
                                        kLHintsExpandY, 3, 3, 3, 3));
   MapSubwindows();
   Resize();
   CenterOnParent();
   Popup();
}

//____________________________________________________________________________
TProofProgressLog::~TProofProgressLog()
{
   // Destructor

   // Cleanup the log object
   SafeDelete(fProofLog);

   // Detach from owner dialog
   if (fDialog) {
      fDialog->fLogWindow = 0;
      fDialog->fProof->Disconnect("LogMessage(const char*,Bool_t)", this,
                                 "LogMessage(const char*,Bool_t)");
   }
}

//____________________________________________________________________________
void TProofProgressLog::Popup()
{
   // Show log window.

   MapWindow();
}

//____________________________________________________________________________
void TProofProgressLog::Clear(Option_t *)
{
   // Clear log window.

   if (fText)
      fText->Clear();
}

//____________________________________________________________________________
void TProofProgressLog::LoadBuffer(const char *buffer)
{
   // Load a text buffer in the window.

   if (fText)
      fText->LoadBuffer(buffer);
}

//____________________________________________________________________________
void TProofProgressLog::LoadFile(const char *file)
{
   // Load a file in the window.

   if (fText)
      fText->LoadFile(file);
}

//____________________________________________________________________________
void TProofProgressLog::AddBuffer(const  char *buffer)
{
   // Add text to the window.

   if (fText) {
      TGText txt;
      txt.LoadBuffer(buffer);
      fText->AddText(&txt);
   }
}

//____________________________________________________________________________
void TProofProgressLog::CloseWindow()
{
   // Handle close button or when closed via window manager action.

   DeleteWindow();
}

//______________________________________________________________________________
void TProofProgressLog::BuildLogList(Bool_t create)
{
   // Build the list of workers. For this, extract the logs and take the names
   // of TProofLogElements

   // Set title
   TString title;
   title.Form("PROOF - Processing logs for session 'undefined'");
   SetWindowName(title.Data());
   SetIconName(title.Data());

   // Create the list-box now
   if (create) {
      if (fLogList) delete fLogList;
      fLogList = new TGListBox(fVworkers);
   } else {
      // Reset
      Int_t nent = fLogList->GetNumberOfEntries();
      fLogList->RemoveEntries(0,nent);
      fLogList->Layout();
   }

   if (fSessionUrl.IsNull()) {
      if (gDebug > 0)
         Info("BuildLogList", "sesssion URL undefined - do nothing");
      return;
   }
   TProofMgr *mgr = TProof::Mgr(fSessionUrl.Data());
   if (!mgr || !mgr->IsValid()) {
      Warning("BuildLogList", "unable open a manager connection to %s",
                              fSessionUrl.Data());
      return;
   }
   if (!(fProofLog = mgr->GetSessionLogs(fSessionIdx,"NR"))) {
      Warning("BuildLogList", "unable to get logs from %s",
                              fSessionUrl.Data());
      return;
   }
   // Set title
   title.Form("PROOF - Processing logs for session '%s', started on %s at %s",
              fProofLog->GetName(), fProofLog->StartTime().AsString(),
              fProofLog->GetTitle());
   SetWindowName(title.Data());
   SetIconName(title.Data());

   TList *elem = fProofLog->GetListOfLogs();
   TIter next(elem);
   TProofLogElem *pe = 0;

   Int_t is = 0;
   TGLBEntry *ent = 0;
   TString buf;
   while ((pe=(TProofLogElem*)next())){
      TUrl url(pe->GetTitle());
      buf.Form("%s %s", pe->GetName(), url.GetHost());
      fLogList->AddEntry(buf.Data(), is);
      if ((ent = fLogList->FindEntry(buf.Data()))) {
         ent->ResetBit(kLogElemFilled);
         ent->ResetBit(kDefaultActive);
         if (!(pe->IsWorker())) ent->SetBit(kDefaultActive);
      }
      is++;
   }

   // Done
   return;
}

//______________________________________________________________________________
void TProofProgressLog::DoLog(Bool_t grep)
{
   // Display the logs

   Clear();

   if (!fGrepText) {
      Warning("DoLog", "no text: do nothing!");
      return;
   }

   TString greptext = fGrepText->GetText();
   Int_t from, to;
   if (fAllLines->IsOn()){
      from = 0;
      to = -1;
   } else {
      from = fLinesFrom->GetIntNumber();
      to = fLinesTo->GetIntNumber();
   }

   // Create the TProofLog instance
   if (!fProofLog) {
      TProofMgr *mgr = 0;
      if ((mgr = TProof::Mgr(fSessionUrl.Data()))) {
         if (!(fProofLog = mgr->GetSessionLogs(fSessionIdx, "NR"))) {
            Warning("DoLog", "unable to instantiate TProofLog for %s",
                             fSessionUrl.Data());
         }
      } else {
         Warning("DoLog", "unable to instantiate a TProofMgr for %s",
                          fSessionUrl.Data());
      }
   }

   // Default is not retrieving
   Bool_t retrieve = kFALSE;
   if (!grep) {
      if (!fFullText ||
          ((fTextType != kRaw && fRawLines->IsOn())   ||
           (fTextType != kStd && !fRawLines->IsOn())) ||
          (fDialog && fDialog->fStatus==TProofProgressDialog::kRunning)) {
         retrieve = kTRUE;
         if (fRawLines->IsOn()) {
            fTextType = kRaw;
         } else {
            fTextType = kStd;
         }
         if (fDialog && fDialog->fStatus != TProofProgressDialog::kRunning)
            fFullText = kTRUE;
      }
   } else {
      retrieve = kTRUE;
      fTextType = kGrep;
      if (fDialog && fDialog->fStatus != TProofProgressDialog::kRunning)
         fFullText = kTRUE;
   }

   // Display now
   if (fProofLog) {
      TList *selected = new TList;
      fLogList->GetSelectedEntries(selected);
      TIter next(selected);
      TGTextLBEntry *selentry;
      Bool_t logonly = fProofLog->LogToBox();
      fProofLog->SetLogToBox(kTRUE);

      fProofLog->Connect("Prt(const char*)", "TProofProgressLog",
                           this, "LogMessage(const char*, Bool_t)");
      while ((selentry=(TGTextLBEntry*)next())){
         TString ord = selentry->GetText()->GetString();
         Int_t is = ord.Index(" ");
         if (is != kNPOS) ord.Remove(is);
         if (retrieve || !selentry->TestBit(kLogElemFilled)) {
            if (fTextType == kGrep) {
               fProofLog->Retrieve(ord.Data(), TProofLog::kGrep, 0, greptext.Data());
            } else if (fTextType == kRaw) {
               fProofLog->Retrieve(ord.Data(), TProofLog::kTrailing, 0, 0);
            } else {
               fProofLog->Retrieve(ord.Data(), TProofLog::kGrep, 0, "-v \"| SvcMsg\"");
            }
            selentry->SetBit(kLogElemFilled);
         }
         fProofLog->Display(ord.Data(), from, to);
      }
      fProofLog->SetLogToBox(logonly);
      fProofLog->Disconnect("Prt(const char*)", this, "LogMessage(const char*, Bool_t)");
      delete selected;
   }
}

//______________________________________________________________________________
void TProofProgressLog::LogMessage(const char *msg, Bool_t all)
{
   // Load/append a log msg in the log frame, if open

   if (all) {
      // load buffer
      LoadBuffer(msg);
   } else {
      // append
      AddBuffer(msg);
   }
}

//______________________________________________________________________________
void TProofProgressLog::SaveToFile()
{
   //Save the logs to a file 
   //Only the name of the file is taken, no expansion

   if (!fProofLog) DoLog();

   // File name: the default is <session-tag>.log
   TString filename = fFileName->GetText();
   if (filename.IsNull() || filename == "<session-tag>.log") {
      filename = (fDialog && fDialog->fProof) ? 
                  TString::Format("%s.log", fDialog->fProof->GetName()) :
                  TString("proof.log");
   }

   TList *selected = new TList;
   fLogList->GetSelectedEntries(selected);
   TIter next(selected);
   TGTextLBEntry *selentry;
   Bool_t writemode=kTRUE;
   const char *option;
   TString ord;
   while ((selentry=(TGTextLBEntry*)next())){
      ord = selentry->GetText()->GetString();
      Int_t isp = ord.Index(' ');
      if (isp != kNPOS) ord.Remove(isp);
      //open the file in "w" mode for the first time
      option = writemode ? "w" : "a";
      fProofLog->Save(ord.Data(), filename.Data(), option);
      writemode=kFALSE;
   }

   Info("SaveToFile", "logs saved to file %s", filename.Data());
   return;
}

//______________________________________________________________________________
void TProofProgressLog::NoLineEntry()
{
   //Enable/disable the line number entry

   if (fAllLines->IsOn()){
      //disable the line number entry
      fLinesFrom->SetState(kFALSE);
      fLinesTo->SetState(kFALSE);
   } else {
      fLinesFrom->SetState(kTRUE);
      fLinesTo->SetState(kTRUE);
   }
}

//______________________________________________________________________________
void TProofProgressLog::Select(Int_t id, Bool_t all)
{
   //actions of select all/clear all button

   Int_t nen = fLogList->GetNumberOfEntries();
   Bool_t sel = id ? 0 : 1;

   TGLBEntry *ent = 0;
   for (Int_t ie=0; ie<nen; ie++) {
      if (all) {
         fLogList->Select(ie, sel);
      } else {
         if ((ent = fLogList->GetEntry(ie))) {
            if (ent->TestBit(kDefaultActive)) fLogList->Select(ie, sel);
         }
      }
   }
}


//______________________________________________________________________________
void TProofProgressLog::Rebuild()
{
   // Rebuild the log info for a new entered session

   // Check if we need to remake the TProofLog object
   Bool_t sameurl = kFALSE;
   TUrl url(fUrlText->GetText());
   TUrl urlref(fSessionUrl.Data());
   if (!strcmp(url.GetHostFQDN(), urlref.GetHostFQDN())) {
      if (url.GetPort() == urlref.GetPort()) {
         if (!strcmp(url.GetUser(), urlref.GetUser())) {
            sameurl = kTRUE;
         }
      }
   }
   Int_t idx = 0;
   if (sameurl) {
      idx = fSessNum->GetIntNumber();
      if (idx == fSessionIdx) {
         Info("Rebuild", "same paremeters {%s, %s}, {%d, %d}: no need to rebuild TProofLog",
                         url.GetUrl(), urlref.GetUrl(), idx, fSessionIdx);
         return;
      }
   }
   // Cleanup current TProofLog
   if (fProofLog) delete fProofLog;

   // Set new parameters
   fSessionUrl = fUrlText->GetText();
   fSessionIdx = idx;

   // Rebuild the list now
   BuildLogList(kFALSE);

   // Select the defaut actives to start with
   Select(0, kFALSE);
   // Redraw
   fLogList->Layout();

   // Done
   return;
}