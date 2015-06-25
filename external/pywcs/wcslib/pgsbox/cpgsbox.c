/*============================================================================

  PGSBOX 4.8 - draw curvilinear coordinate axes for PGPLOT.
  Copyright (C) 1997-2011, Mark Calabretta

  This file is part of PGSBOX.

  PGSBOX is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  PGSBOX is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with PGSBOX.  If not, see <http://www.gnu.org/licenses/>.

  Correspondence concerning PGSBOX may be directed to:
    Internet email: mcalabre@atnf.csiro.au
    Postal address: Dr. Mark Calabretta
                    Australia Telescope National Facility, CSIRO
                    PO Box 76
                    Epping NSW 1710
                    AUSTRALIA

  Author: Mark Calabretta, Australia Telescope National Facility
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: cpgsbox.c,v 1.1.2.1 2013/12/14 18:47:06 jasercio Exp $
*===========================================================================*/

#include <string.h>
#include "cpgsbox.h"

/* Fortran name mangling. */
#include <wcsconfig_f77.h>
#define pgsbok_ F77_FUNC(pgsbok, PGSBOK)
#define pglbok_ F77_FUNC(pglbok, PGLBOK)

void pgsbok_(const float blc[2], const float trc[2], char idents[3][80],
             const char opt[2], const int *labctl, const int *labden,
             const int ci[7], const int gcode[2], const double *tiklen,
             const int *ng1, const double *grid1, const int *ng2,
             const double *grid2, const int *doeq, nlfunc_t nlfunc,
             const int *nlc, const int *nli, const int *nld, char *nlcprm,
             int *nliprm, double *nldprm, const int *nc, int *ic,
             double cache[][4], int *ierr);

void pglbok_(char idents[3][80], const char opt[2], const int *labctl,
             const int *labden, const int ci[7], const int gcode[2],
             const double *tiklen, const int *ng1, const double *grid1,
             const int *ng2, const double *grid2, const int *doeq,
             const int *nc, int *ic, double cache[][4], int *ierr);

void cpgsbox(
  const float blc[2],
  const float trc[2],
  char (*idents)[80],
  const char opt[2],
  int labctl,
  int labden,
  const int ci[7],
  const int gcode[2],
  double tiklen,
  int ng1,
  const double *grid1,
  int ng2,
  const double *grid2,
  int doeq,
  nlfunc_t nlfunc,
  int nlc,
  int nli,
  int nld,
  char nlcprm[],
  int nliprm[],
  double nldprm[],
  int nc,
  int *ic,
  double cache[][4],
  int *ierr)

{
  char ids[3][80];
  int j, k;

  /* Convert variable length strings to fixed-length char arrays. */
  k = 0;
  for (j = 0; j < 3; j++) {
    if (strlen(idents[j]) > 80) {
      strncpy(ids[j], idents[j], 80);
    } else {
      strcpy(ids[j], idents[j]);
      for (k = strlen(idents[j]); k < 80; k++) {
        ids[j][k] = ' ';
      }
    }
  }

  pgsbok_(blc, trc, ids, opt, &labctl, &labden, ci, gcode, &tiklen, &ng1,
          grid1, &ng2, grid2, &doeq, nlfunc, &nlc, &nli, &nld, nlcprm, nliprm,
          nldprm, &nc, ic, cache, ierr);
  return;
}

/*==========================================================================*/

void cpglbox(
  char (*idents)[80],
  const char opt[2],
  int labctl,
  int labden,
  const int ci[7],
  const int gcode[2],
  double tiklen,
  int ng1,
  const double *grid1,
  int ng2,
  const double *grid2,
  int doeq,
  int nc,
  int *ic,
  double cache[][4],
  int *ierr)

{
  char ids[3][80];
  int j, k;

  /* Convert variable length strings to fixed-length char arrays. */
  k = 0;
  for (j = 0; j < 3; j++) {
    if (strlen(idents[j]) > 80) {
      strncpy(ids[j], idents[j], 80);
    } else {
      strcpy(ids[j], idents[j]);
      for (k = strlen(idents[j]); k < 80; k++) {
        ids[j][k] = ' ';
      }
    }
  }

  pglbok_(ids, opt, &labctl, &labden, ci, gcode, &tiklen, &ng1, grid1, &ng2,
          grid2, &doeq, &nc, ic, cache, ierr);
  return;
}