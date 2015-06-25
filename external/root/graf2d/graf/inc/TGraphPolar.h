// @(#)root/graf:$Id: TGraphPolar.h,v 1.1.1.2 2010/04/07 18:31:50 elwinter Exp $
// Author: Sebastian Boser, 02/02/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphPolar
#define ROOT_TGraphPolar

#ifndef ROOT_TGraphErrors
#include "TGraphErrors.h"
#endif
#ifndef ROOT_Riosfwd
#include "Riosfwd.h"
#endif
#ifndef ROOT_TAttText
#include "TAttText.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif

#include "TGraphPolargram.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraphPolar                                                          //
//                                                                      //
// Polar graph graphics class.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TGraphPolar: public TGraphErrors {

private:
   Bool_t fOptionAxis;          // Force drawing of new coord system

protected:
   TGraphPolargram* fPolargram; // The polar coordinates system
   Double_t* fXpol;             // [fNpoints] points in polar coordinates
   Double_t* fYpol;             // [fNpoints] points in polar coordinates


public:
   TGraphPolar();
   TGraphPolar(Int_t n, const Double_t* x=0, const Double_t* y=0,
                        const Double_t* ex=0, const Double_t* ey=0);
   virtual ~TGraphPolar();

   TGraphPolargram *GetPolargram() {return fPolargram;};

   void             Draw(Option_t* options = "");
   Bool_t           GetOptionAxis() {return fOptionAxis;};
   void             SetMaxRadial(Double_t maximum = 1); //*MENU*
   void             SetMinRadial(Double_t minimum = 0); //*MENU*
   void             SetMaximum(Double_t maximum = 1) {SetMaxRadial(maximum);}
   void             SetMinimum(Double_t minimum = 0) {SetMinRadial(minimum);}
   void             SetMaxPolar(Double_t maximum = 6.28318530717958623); //*MENU*
   void             SetMinPolar(Double_t minimum = 0); //*MENU*
   void             SetOptionAxis(Bool_t opt) {fOptionAxis = opt;};
   void             SetPolargram(TGraphPolargram *p) {fPolargram = p;};
   Double_t        *GetXpol();
   Double_t        *GetYpol();

   ClassDef(TGraphPolar,1); // Polar graph
};

#endif