*=======================================================================
*
* WCSLIB 4.8 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2011, Mark Calabretta
*
* This file is part of WCSLIB.
*
* WCSLIB is free software: you can redistribute it and/or modify it
* under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* WCSLIB is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
* License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with WCSLIB.  If not, see http://www.gnu.org/licenses.
*
* Correspondence concerning WCSLIB may be directed to:
*   Internet email: mcalabre@atnf.csiro.au
*   Postal address: Dr. Mark Calabretta
*                   Australia Telescope National Facility, CSIRO
*                   PO Box 76
*                   Epping NSW 1710
*                   AUSTRALIA
*
* Author: Mark Calabretta, Australia Telescope National Facility
* http://www.atnf.csiro.au/~mcalabre/index.html
* $Id: prj.inc,v 4.8.1.1 2011/08/15 08:07:07 cal103 Exp cal103 $
*=======================================================================

*     Functions.
      EXTERNAL  PRJFREE, PRJGET, PRJGTC, PRJGTD, PRJGTI, PRJINI, PRJPRT,
     :          PRJPTC,  PRJPTD, PRJPTI, PRJPUT, PRJS2X, PRJSET, PRJX2S
      INTEGER   PRJFREE, PRJGET, PRJGTC, PRJGTD, PRJGTI, PRJINI, PRJPRT,
     :          PRJPTC,  PRJPTD, PRJPTI, PRJPUT, PRJS2X, PRJSET, PRJX2S

      EXTERNAL  AZPSET, AZPX2S, AZPS2X,    SFLSET, SFLX2S, SFLS2X,
     :          SZPSET, SZPX2S, SZPS2X,    PARSET, PARX2S, PARS2X,
     :          TANSET, TANX2S, TANS2X,    MOLSET, MOLX2S, MOLS2X,
     :          STGSET, STGX2S, STGS2X,    AITSET, AITX2S, AITS2X,
     :          SINSET, SINX2S, SINS2X,    COPSET, COPX2S, COPS2X,
     :          ARCSET, ARCX2S, ARCS2X,    COESET, COEX2S, COES2X,
     :          ZPNSET, ZPNX2S, ZPNS2X,    CODSET, CODX2S, CODS2X,
     :          ZEASET, ZEAX2S, ZEAS2X,    COOSET, COOX2S, COOS2X,
     :          AIRSET, AIRX2S, AIRS2X,    BONSET, BONX2S, BONS2X,
     :          CYPSET, CYPX2S, CYPS2X,    PCOSET, PCOX2S, PCOS2X,
     :          CEASET, CEAX2S, CEAS2X,    TSCSET, TSCX2S, TSCS2X,
     :          CARSET, CARX2S, CARS2X,    CSCSET, CSCX2S, CSCS2X,
     :          MERSET, MERX2S, MERS2X,    QSCSET, QSCX2S, QSCS2X
      INTEGER   AZPSET, AZPX2S, AZPS2X,    SFLSET, SFLX2S, SFLS2X,
     :          SZPSET, SZPX2S, SZPS2X,    PARSET, PARX2S, PARS2X,
     :          TANSET, TANX2S, TANS2X,    MOLSET, MOLX2S, MOLS2X,
     :          STGSET, STGX2S, STGS2X,    AITSET, AITX2S, AITS2X,
     :          SINSET, SINX2S, SINS2X,    COPSET, COPX2S, COPS2X,
     :          ARCSET, ARCX2S, ARCS2X,    COESET, COEX2S, COES2X,
     :          ZPNSET, ZPNX2S, ZPNS2X,    CODSET, CODX2S, CODS2X,
     :          ZEASET, ZEAX2S, ZEAS2X,    COOSET, COOX2S, COOS2X,
     :          AIRSET, AIRX2S, AIRS2X,    BONSET, BONX2S, BONS2X,
     :          CYPSET, CYPX2S, CYPS2X,    PCOSET, PCOX2S, PCOS2X,
     :          CEASET, CEAX2S, CEAS2X,    TSCSET, TSCX2S, TSCS2X,
     :          CARSET, CARX2S, CARS2X,    CSCSET, CSCX2S, CSCS2X,
     :          MERSET, MERX2S, MERS2X,    QSCSET, QSCX2S, QSCS2X

*     Length of the PRJPRM data structure (INTEGER array) on 64-bit
*     machines.  Only needs to be 116 on 32-bit machines.
      INTEGER   PRJLEN
      PARAMETER (PRJLEN = 120)

*     Number of projection parameters supported by WCSLIB, 0 to PVN-1.
      INTEGER   PRJ_PVN
      PARAMETER (PRJ_PVN = 30)

*     Codes for PRJ data structure elements used by PRJPUT and PRJGET.
      INTEGER   PRJ_BOUNDS, PRJ_CODE, PRJ_FLAG, PRJ_PHI0, PRJ_PV,
     :          PRJ_R0, PRJ_THETA0

      PARAMETER (PRJ_FLAG      = 100)
      PARAMETER (PRJ_CODE      = 101)
      PARAMETER (PRJ_R0        = 102)
      PARAMETER (PRJ_PV        = 103)
      PARAMETER (PRJ_PHI0      = 104)
      PARAMETER (PRJ_THETA0    = 105)
      PARAMETER (PRJ_BOUNDS    = 106)

*     Codes for PRJ data structure elements used by PRJGET (only).
      INTEGER   PRJ_CATEGORY, PRJ_CONFORMAL, PRJ_ERR, PRJ_GLOBAL,
     :          PRJ_DIVERGENT, PRJ_EQUIAREAL, PRJ_N, PRJ_NAME,
     :          PRJ_PVRANGE, PRJ_SIMPLEZEN, PRJ_W, PRJ_X0, PRJ_Y0

      PARAMETER (PRJ_NAME      = 200)
      PARAMETER (PRJ_CATEGORY  = 201)
      PARAMETER (PRJ_PVRANGE   = 202)
      PARAMETER (PRJ_SIMPLEZEN = 203)
      PARAMETER (PRJ_EQUIAREAL = 204)
      PARAMETER (PRJ_CONFORMAL = 205)
      PARAMETER (PRJ_GLOBAL    = 206)
      PARAMETER (PRJ_DIVERGENT = 207)
      PARAMETER (PRJ_X0        = 208)
      PARAMETER (PRJ_Y0        = 209)
      PARAMETER (PRJ_ERR       = 210)
      PARAMETER (PRJ_W         = 211)
      PARAMETER (PRJ_N         = 212)

*     Projection categories.
      INTEGER   PRJ_CONIC, PRJ_CONVENTIONAL, PRJ_CYLINDRICAL,
     :          PRJ_HEALPIX, PRJ_POLYCONIC, PRJ_PSEUDOCYLINDRICAL,
     :          PRJ_QUADCUBE, PRJ_ZENITHAL

      PARAMETER (PRJ_ZENITHAL          = 1)
      PARAMETER (PRJ_CYLINDRICAL       = 2)
      PARAMETER (PRJ_PSEUDOCYLINDRICAL = 3)
      PARAMETER (PRJ_CONVENTIONAL      = 4)
      PARAMETER (PRJ_CONIC             = 5)
      PARAMETER (PRJ_POLYCONIC         = 6)
      PARAMETER (PRJ_QUADCUBE          = 7)
      PARAMETER (PRJ_HEALPIX           = 8)

*     Error codes and messages.
      INTEGER   PRJERR_BAD_PARAM, PRJERR_BAD_PIX, PRJERR_BAD_WORLD,
     :          PRJERR_NULL_POINTER, PRJERR_SUCCESS

      PARAMETER (PRJERR_SUCCESS      = 0)
      PARAMETER (PRJERR_NULL_POINTER = 1)
      PARAMETER (PRJERR_BAD_PARAM    = 2)
      PARAMETER (PRJERR_BAD_PIX      = 3)
      PARAMETER (PRJERR_BAD_WORLD    = 4)

      CHARACTER PRJ_ERRMSG(0:4)*80
      COMMON /PRJ_DATA/ PRJ_ERRMSG
