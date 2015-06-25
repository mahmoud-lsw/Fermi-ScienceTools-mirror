// @(#)root/ged:$Id: TAttLineEditor.h,v 1.1.1.1 2009/04/05 20:43:50 elwinter Exp $
// Author: Ilka  Antcheva 10/05/04

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAttLineEditor
#define ROOT_TAttLineEditor

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TAttLineEditor                                                      //
//                                                                      //
//  Implements GUI for editing line attributes.                         //                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGedFrame
#include "TGedFrame.h"
#endif

class TGLineStyleComboBox;
class TGLineWidthComboBox;
class TGColorSelect;
class TAttLine;

class TAttLineEditor : public TGedFrame {

protected:
   TAttLine             *fAttLine;          // line attribute object
   TGLineStyleComboBox  *fStyleCombo;       // line style combo box
   TGLineWidthComboBox  *fWidthCombo;       // line width combo box
   TGColorSelect        *fColorSelect;      // line color widget

   virtual void   ConnectSignals2Slots();

public:
   TAttLineEditor(const TGWindow *p = 0, 
                  Int_t width = 140, Int_t height = 30,
                  UInt_t options = kChildFrame,
                  Pixel_t back = GetDefaultFrameBackground());
   virtual ~TAttLineEditor();

   virtual void   SetModel(TObject* obj);
   virtual void   DoLineColor(Pixel_t color);
   virtual void   DoLineStyle(Int_t style);
   virtual void   DoLineWidth(Int_t width);

   ClassDef(TAttLineEditor,0)  // GUI for editing line attributes
};

#endif
