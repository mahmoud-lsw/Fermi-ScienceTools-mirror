// @(#)root/eve:$Id: TEveTriangleSetGL.cxx,v 1.1.1.2 2013/01/23 16:01:40 areustle Exp $
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveTriangleSetGL.h"
#include "TEveTriangleSet.h"
#include "TVector3.h"

#include "TGLIncludes.h"

//==============================================================================
// TEveTriangleSetGL
//==============================================================================

//______________________________________________________________________________
//
// GL-renderer for TEveTriangleSet class.
//
// See also: TGLObject, TGLLogicalShape.

ClassImp(TEveTriangleSetGL);

//______________________________________________________________________________
TEveTriangleSetGL::TEveTriangleSetGL() : TGLObject(), fM(0)
{
   // Constructor.

   // fDLCache = false; // Disable display list.
   fMultiColor = kTRUE;
}

//______________________________________________________________________________
TEveTriangleSetGL::~TEveTriangleSetGL()
{
   // Destructor.
}

/******************************************************************************/

//______________________________________________________________________________
Bool_t TEveTriangleSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
   // Set model object.

   fM = SetModelDynCast<TEveTriangleSet>(obj);
   return kTRUE;
}

//______________________________________________________________________________
void TEveTriangleSetGL::SetBBox()
{
   // Set bounding-box from the model.

   // !! This ok if master sub-classed from TAttBBox
   SetAxisAlignedBBox(((TEveTriangleSet*)fExternalObj)->AssertBBox());
}

/******************************************************************************/

//______________________________________________________________________________
void TEveTriangleSetGL::DirectDraw(TGLRnrCtx & /*rnrCtx*/) const
{
   // Low-level GL rendering.

   TEveTriangleSet& refTS = *fM;
   Bool_t isScaled = refTS.RefMainTrans().IsScale();

   GLint ex_shade_model;
   glGetIntegerv(GL_SHADE_MODEL, &ex_shade_model);
   glShadeModel(GL_FLAT);

   glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

   glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);
   glDisable(GL_CULL_FACE);
   if (isScaled) glEnable(GL_NORMALIZE);
   glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
   glVertexPointer(3, GL_FLOAT, 0, refTS.fVerts);
   glEnableClientState(GL_VERTEX_ARRAY);

   Int_t*   tng = refTS.fTrings;
   Float_t* nrm = refTS.fTringNorms;
   UChar_t* col = refTS.fTringCols;

   TVector3 e1, e2, n;

   glBegin(GL_TRIANGLES);
   for(Int_t t=0; t<refTS.fNTrings; ++t) {
      if (nrm) {
         glNormal3fv(nrm); nrm += 3;
      } else {
         Float_t* v0 = refTS.Vertex(tng[0]);
         Float_t* v1 = refTS.Vertex(tng[1]);
         Float_t* v2 = refTS.Vertex(tng[2]);
         e1.SetXYZ(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
         e2.SetXYZ(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
         n = e1.Cross(e2);
         if (!isScaled) n.SetMag(1);
         glNormal3d(n.x(), n.y(), n.z());
      }
      if (col) {
         TGLUtil::Color3ubv(col);  col += 3;
      }
      glArrayElement(tng[0]);
      glArrayElement(tng[1]);
      glArrayElement(tng[2]);
      tng += 3;
   }
   glEnd();

   glPopClientAttrib();
   glPopAttrib();
   glShadeModel(ex_shade_model);
}