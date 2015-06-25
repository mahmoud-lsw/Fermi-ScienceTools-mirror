/* drmnfb_routines.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c/f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b9 = 0.;
static integer c__3 = 3;
static doublereal c_b25 = .33333333333333331;
static doublereal c_b26 = -.66666666666666663;
static doublereal c_b28 = .66666666666666663;

/* Subroutine */ int drmnfb_(doublereal *b, doublereal *d__, doublereal *fx, 
	integer *iv, integer *liv, integer *lv, integer *p, doublereal *v, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, w, g1, iv1, ipi, alpha, alpha0;
    extern /* Subroutine */ int ds3grd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dv7scp_(integer *, doublereal *, doublereal *)
	    , drmngb_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), divset_(integer *, integer *, integer *, integer *,
	     doublereal *);


/*  ***  ITERATION DRIVER FOR  DMNF... */
/*  ***  MINIMIZE GENERAL UNCONSTRAINED OBJECTIVE FUNCTION USING */
/*  ***  FINITE-DIFFERENCE GRADIENTS AND SECANT HESSIAN APPROXIMATIONS. */

/*     DIMENSION IV(59 + P), V(77 + P*(P+23)/2) */

/*  ***  PURPOSE  *** */

/*        THIS ROUTINE INTERACTS WITH SUBROUTINE  DRMNGB  IN AN ATTEMPT */
/*     TO FIND AN P-VECTOR  X*  THAT MINIMIZES THE (UNCONSTRAINED) */
/*     OBJECTIVE FUNCTION  FX = F(X)  COMPUTED BY THE CALLER.  (OFTEN */
/*     THE  X*  FOUND IS A LOCAL MINIMIZER RATHER THAN A GLOBAL ONE.) */

/*  ***  PARAMETERS  *** */

/*        THE PARAMETERS FOR DRMNFB ARE THE SAME AS THOSE FOR  DMNG */
/*     (WHICH SEE), EXCEPT THAT CALCF, CALCG, UIPARM, URPARM, AND UFPARM */
/*     ARE OMITTED, AND A PARAMETER  FX  FOR THE OBJECTIVE FUNCTION */
/*     VALUE AT X IS ADDED.  INSTEAD OF CALLING CALCG TO OBTAIN THE */
/*     GRADIENT OF THE OBJECTIVE FUNCTION AT X, DRMNFB CALLS DS3GRD, */
/*     WHICH COMPUTES AN APPROXIMATION TO THE GRADIENT BY FINITE */
/*     (FORWARD AND CENTRAL) DIFFERENCES USING THE METHOD OF REF. 1. */
/*     THE FOLLOWING INPUT COMPONENT IS OF INTEREST IN THIS REGARD */
/*     (AND IS NOT DESCRIBED IN  DMNG). */

/* V(ETA0)..... V(42) IS AN ESTIMATED BOUND ON THE RELATIVE ERROR IN THE */
/*             OBJECTIVE FUNCTION VALUE COMPUTED BY CALCF... */
/*                  (TRUE VALUE) = (COMPUTED VALUE) * (1 + E), */
/*             WHERE ABS(E) .LE. V(ETA0).  DEFAULT = MACHEP * 10**3, */
/*             WHERE MACHEP IS THE UNIT ROUNDOFF. */

/*        THE OUTPUT VALUES IV(NFCALL) AND IV(NGCALL) HAVE DIFFERENT */
/*     MEANINGS FOR  DMNF THAN FOR  DMNG... */

/* IV(NFCALL)... IV(6) IS THE NUMBER OF CALLS SO FAR MADE ON CALCF (I.E., */
/*             FUNCTION EVALUATIONS) EXCLUDING THOSE MADE ONLY FOR */
/*             COMPUTING GRADIENTS.  THE INPUT VALUE IV(MXFCAL) IS A */
/*             LIMIT ON IV(NFCALL). */
/* IV(NGCALL)... IV(30) IS THE NUMBER OF FUNCTION EVALUATIONS MADE ONLY */
/*             FOR COMPUTING GRADIENTS.  THE TOTAL NUMBER OF FUNCTION */
/*             EVALUATIONS IS THUS  IV(NFCALL) + IV(NGCALL). */

/*  ***  REFERENCES  *** */

/* 1. STEWART, G.W. (1967), A MODIFICATION OF DAVIDON*S MINIMIZATION */
/*        METHOD TO ACCEPT DIFFERENCE APPROXIMATIONS OF DERIVATIVES, */
/*        J. ASSOC. COMPUT. MACH. 14, PP. 72-83. */
/* . */
/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (AUGUST 1982). */

/* ----------------------------  DECLARATIONS  --------------------------- */


/* DIVSET.... SUPPLIES DEFAULT PARAMETER VALUES. */
/* DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* DS3GRD... COMPUTES FINITE-DIFFERENCE GRADIENT APPROXIMATION. */
/* DRMNGB... REVERSE-COMMUNICATION ROUTINE THAT DOES  DMNGB ALGORITHM. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */


/*  ***  SUBSCRIPTS FOR IV   *** */


/* /6 */
/*     DATA ETA0/42/, F/10/, G/28/, LMAT/42/, NEXTV/47/, NGCALL/30/, */
/*    1     NITER/31/, PERM/58/, SGIRC/57/, TOOBIG/2/, VNEED/4/ */
/* /7 */
/* / */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --d__;
    b -= 3;

    /* Function Body */
    iv1 = iv[1];
    if (iv1 == 1) {
	goto L10;
    }
    if (iv1 == 2) {
	goto L50;
    }
    if (iv[1] == 0) {
	divset_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 == 12 || iv1 == 13) {
	iv[4] = iv[4] + (*p << 1) + 6;
    }
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    g1 = 1;
    if (iv1 == 12) {
	iv[1] = 13;
    }
    goto L20;

L10:
    g1 = iv[28];

L20:
    drmngb_(&b[3], &d__[1], fx, &v[g1], &iv[1], liv, lv, p, &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L999;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L80;
    }

/*  ***  COMPUTE GRADIENT  *** */

L30:
    if (iv[31] == 0) {
	dv7scp_(p, &v[g1], &c_b9);
    }
    j = iv[42];
    alpha0 = g1 - *p - 1;
    ipi = iv[58];
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = alpha0 + iv[ipi];
	v[k] = dd7tpr_(&i__, &v[j], &v[j]);
	++ipi;
	j += i__;
/* L40: */
    }
/*     ***  UNDO INCREMENT OF IV(NGCALL) DONE BY DRMNGB  *** */
    --iv[30];
/*     ***  STORE RETURN CODE FROM DS3GRD IN IV(SGIRC)  *** */
    iv[57] = 0;
/*     ***  X MAY HAVE BEEN RESTORED, SO COPY BACK FX... *** */
    *fx = v[10];
    goto L60;

/*     ***  GRADIENT LOOP  *** */

L50:
    if (iv[2] != 0) {
	goto L10;
    }

L60:
    g1 = iv[28];
    alpha = g1 - *p;
    w = alpha - 6;
    ds3grd_(&v[alpha], &b[3], &d__[1], &v[42], fx, &v[g1], &iv[57], p, &v[w], 
	    &x[1]);
    i__ = iv[57];
    if (i__ == 0) {
	goto L10;
    }
    if (i__ <= *p) {
	goto L70;
    }
    iv[2] = 1;
    goto L10;

L70:
    ++iv[30];
    goto L999;

L80:
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION  *** */

    iv[28] = iv[47] + *p + 6;
    iv[47] = iv[28] + *p;
    if (iv1 != 13) {
	goto L10;
    }

L999:
    return 0;
/*  ***  LAST CARD OF DRMNFB FOLLOWS  *** */
} /* drmnfb_ */

/* Subroutine */ int ds3grd_(doublereal *alpha, doublereal *b, doublereal *
	d__, doublereal *eta0, doublereal *fx, doublereal *g, integer *irc, 
	integer *p, doublereal *w, doublereal *x)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal h0, gi, xi, aai, agi, eta, afx, axi;
    static logical hit;
    static doublereal xih, hmin;
    extern doublereal dr7mdc_(integer *);
    static doublereal machep, alphai, axibar, afxeta, discon;


/*  ***  COMPUTE FINITE DIFFERENCE GRADIENT BY STWEART*S SCHEME  *** */

/*     ***  PARAMETERS  *** */


/* ....................................................................... */

/*     ***  PURPOSE  *** */

/*        THIS SUBROUTINE USES AN EMBELLISHED FORM OF THE FINITE-DIFFER- */
/*     ENCE SCHEME PROPOSED BY STEWART (REF. 1) TO APPROXIMATE THE */
/*     GRADIENT OF THE FUNCTION F(X), WHOSE VALUES ARE SUPPLIED BY */
/*     REVERSE COMMUNICATION. */

/*     ***  PARAMETER DESCRIPTION  *** */

/*  ALPHA IN  (APPROXIMATE) DIAGONAL ELEMENTS OF THE HESSIAN OF F(X). */
/*      B IN  ARRAY OF SIMPLE LOWER AND UPPER BOUNDS ON X.  X MUST */
/*             SATISFY B(1,I) .LE. X(I) .LE. B(2,I), I = 1(1)P. */
/*             FOR ALL I WITH B(1,I) .GE. B(2,I), DS3GRD SIMPLY */
/*             SETS G(I) TO 0. */
/*      D IN  SCALE VECTOR SUCH THAT D(I)*X(I), I = 1,...,P, ARE IN */
/*             COMPARABLE UNITS. */
/*   ETA0 IN  ESTIMATED BOUND ON RELATIVE ERROR IN THE FUNCTION VALUE... */
/*             (TRUE VALUE) = (COMPUTED VALUE)*(1+E),   WHERE */
/*             ABS(E) .LE. ETA0. */
/*     FX I/O ON INPUT,  FX  MUST BE THE COMPUTED VALUE OF F(X).  ON */
/*             OUTPUT WITH IRC = 0, FX HAS BEEN RESTORED TO ITS ORIGINAL */
/*             VALUE, THE ONE IT HAD WHEN DS3GRD WAS LAST CALLED WITH */
/*             IRC = 0. */
/*      G I/O ON INPUT WITH IRC = 0, G SHOULD CONTAIN AN APPROXIMATION */
/*             TO THE GRADIENT OF F NEAR X, E.G., THE GRADIENT AT THE */
/*             PREVIOUS ITERATE.  WHEN DS3GRD RETURNS WITH IRC = 0, G IS */
/*             THE DESIRED FINITE-DIFFERENCE APPROXIMATION TO THE */
/*             GRADIENT AT X. */
/*    IRC I/O INPUT/RETURN CODE... BEFORE THE VERY FIRST CALL ON DS3GRD, */
/*             THE CALLER MUST SET IRC TO 0.  WHENEVER DS3GRD RETURNS A */
/*             NONZERO VALUE (OF AT MOST P) FOR IRC, IT HAS PERTURBED */
/*             SOME COMPONENT OF X... THE CALLER SHOULD EVALUATE F(X) */
/*             AND CALL DS3GRD AGAIN WITH FX = F(X).  IF B PREVENTS */
/*             ESTIMATING G(I) I.E., IF THERE IS AN I WITH */
/*             B(1,I) .LT. B(2,I) BUT WITH B(1,I) SO CLOSE TO B(2,I) */
/*             THAT THE FINITE-DIFFERENCING STEPS CANNOT BE CHOSEN, */
/*             THEN DS3GRD RETURNS WITH IRC .GT. P. */
/*      P IN  THE NUMBER OF VARIABLES (COMPONENTS OF X) ON WHICH F */
/*             DEPENDS. */
/*      X I/O ON INPUT WITH IRC = 0, X IS THE POINT AT WHICH THE */
/*             GRADIENT OF F IS DESIRED.  ON OUTPUT WITH IRC NONZERO, X */
/*             IS THE POINT AT WHICH F SHOULD BE EVALUATED.  ON OUTPUT */
/*             WITH IRC = 0, X HAS BEEN RESTORED TO ITS ORIGINAL VALUE */
/*             (THE ONE IT HAD WHEN DS3GRD WAS LAST CALLED WITH IRC = 0) */
/*             AND G CONTAINS THE DESIRED GRADIENT APPROXIMATION. */
/*      W I/O WORK VECTOR OF LENGTH 6 IN WHICH DS3GRD SAVES CERTAIN */
/*             QUANTITIES WHILE THE CALLER IS EVALUATING F(X) AT A */
/*             PERTURBED X. */

/*     ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*        THIS ROUTINE IS INTENDED FOR USE WITH QUASI-NEWTON ROUTINES */
/*     FOR UNCONSTRAINED MINIMIZATION (IN WHICH CASE  ALPHA  COMES FROM */
/*     THE DIAGONAL OF THE QUASI-NEWTON HESSIAN APPROXIMATION). */

/*     ***  ALGORITHM NOTES  *** */

/*        THIS CODE DEPARTS FROM THE SCHEME PROPOSED BY STEWART (REF. 1) */
/*     IN ITS GUARDING AGAINST OVERLY LARGE OR SMALL STEP SIZES AND ITS */
/*     HANDLING OF SPECIAL CASES (SUCH AS ZERO COMPONENTS OF ALPHA OR G). */

/*     ***  REFERENCES  *** */

/* 1. STEWART, G.W. (1967), A MODIFICATION OF DAVIDON*S MINIMIZATION */
/*        METHOD TO ACCEPT DIFFERENCE APPROXIMATIONS OF DERIVATIVES, */
/*        J. ASSOC. COMPUT. MACH. 14, PP. 72-83. */

/*     ***  HISTORY  *** */

/*     DESIGNED AND CODED BY DAVID M. GAY (SUMMER 1977/SUMMER 1980). */

/*     ***  GENERAL  *** */

/*        THIS ROUTINE WAS PREPARED IN CONNECTION WITH WORK SUPPORTED BY */
/*     THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS76-00324 AND */
/*     MCS-7906671. */

/* ....................................................................... */

/*     *****  EXTERNAL FUNCTION  ***** */

/* DR7MDC... RETURNS MACHINE-DEPENDENT CONSTANTS. */

/*     ***** INTRINSIC FUNCTIONS ***** */
/* /+ */
/* / */
/*     ***** LOCAL VARIABLES ***** */


/* /6 */
/*     DATA C2000/2.0D+3/, FOUR/4.0D+0/, HMAX0/0.02D+0/, HMIN0/5.0D+1/, */
/*    1     ONE/1.0D+0/, P002/0.002D+0/, THREE/3.0D+0/, */
/*    2     TWO/2.0D+0/, ZERO/0.0D+0/ */
/* /7 */
/* / */
/* /6 */
/*     DATA FH/3/, FX0/4/, HSAVE/5/, XISAVE/6/ */
/* /7 */
/* / */

/* ---------------------------------  BODY  ------------------------------ */

    /* Parameter adjustments */
    --x;
    --g;
    --d__;
    b -= 3;
    --alpha;
    --w;

    /* Function Body */
    if (*irc < 0) {
	goto L80;
    } else if (*irc == 0) {
	goto L10;
    } else {
	goto L210;
    }

/*     ***  FRESH START -- GET MACHINE-DEPENDENT CONSTANTS  *** */

/*     STORE MACHEP IN W(1) AND H0 IN W(2), WHERE MACHEP IS THE UNIT */
/*     ROUNDOFF (THE SMALLEST POSITIVE NUMBER SUCH THAT */
/*     1 + MACHEP .GT. 1  AND  1 - MACHEP .LT. 1),  AND  H0 IS THE */
/*     SQUARE-ROOT OF MACHEP. */

L10:
    w[1] = dr7mdc_(&c__3);
    w[2] = sqrt(w[1]);

    w[4] = *fx;

/*     ***  INCREMENT  I  AND START COMPUTING  G(I)  *** */

L20:
    i__ = abs(*irc) + 1;
    if (i__ > *p) {
	goto L220;
    }
    *irc = i__;
    if (b[(i__ << 1) + 1] < b[(i__ << 1) + 2]) {
	goto L30;
    }
    g[i__] = 0.;
    goto L20;
L30:
    afx = abs(w[4]);
    machep = w[1];
    h0 = w[2];
    hmin = machep * 50.;
    xi = x[i__];
    w[6] = xi;
    axi = abs(xi);
/* Computing MAX */
    d__1 = axi, d__2 = 1. / d__[i__];
    axibar = max(d__1,d__2);
    gi = g[i__];
    agi = abs(gi);
    eta = abs(*eta0);
    if (afx > 0.) {
/* Computing MAX */
	d__1 = eta, d__2 = agi * axi * machep / afx;
	eta = max(d__1,d__2);
    }
    alphai = alpha[i__];
    if (alphai == 0.) {
	goto L130;
    }
    if (gi == 0. || *fx == 0.) {
	goto L140;
    }
    afxeta = afx * eta;
    aai = abs(alphai);

/*        *** COMPUTE H = STEWART*S FORWARD-DIFFERENCE STEP SIZE. */

/* Computing 2nd power */
    d__1 = gi;
    if (d__1 * d__1 <= afxeta * aai) {
	goto L40;
    }
    h__ = sqrt(afxeta / aai) * 2.;
    h__ *= 1. - aai * h__ / (aai * 3. * h__ + agi * 4.);
    goto L50;
/* 40      H = TWO*(AFXETA*AGI/(AAI**2))**(ONE/THREE) */
L40:
    d__1 = afxeta * agi;
    h__ = pow_dd(&d__1, &c_b25) * 2. * pow_dd(&aai, &c_b26);
    h__ *= 1. - agi * 2. / (aai * 3. * h__ + agi * 4.);

/*        ***  ENSURE THAT  H  IS NOT INSIGNIFICANTLY SMALL  *** */

L50:
/* Computing MAX */
    d__1 = h__, d__2 = hmin * axibar;
    h__ = max(d__1,d__2);

/*        *** USE FORWARD DIFFERENCE IF BOUND ON TRUNCATION ERROR IS AT */
/*        *** MOST 10**-3. */

    if (aai * h__ <= agi * .002) {
	goto L120;
    }

/*        *** COMPUTE H = STEWART*S STEP FOR CENTRAL DIFFERENCE. */

    discon = afxeta * 2e3;
/* Computing 2nd power */
    d__1 = gi;
    h__ = discon / (agi + sqrt(d__1 * d__1 + aai * discon));

/*        ***  ENSURE THAT  H  IS NEITHER TOO SMALL NOR TOO BIG  *** */

/* Computing MAX */
    d__1 = h__, d__2 = hmin * axibar;
    h__ = max(d__1,d__2);
    if (h__ >= axibar * .02) {
	h__ = axibar * pow_dd(&h0, &c_b28);
    }

/*        ***  COMPUTE CENTRAL DIFFERENCE  *** */

    xih = xi + h__;
    if (xi - h__ < b[(i__ << 1) + 1]) {
	goto L60;
    }
    *irc = -i__;
    if (xih <= b[(i__ << 1) + 2]) {
	goto L200;
    }
    h__ = -h__;
    xih = xi + h__;
    if (xi + h__ * 2. < b[(i__ << 1) + 1]) {
	goto L190;
    }
    goto L70;
L60:
    if (xi + h__ * 2. > b[(i__ << 1) + 2]) {
	goto L190;
    }
/*        *** MUST DO OFF-SIDE CENTRAL DIFFERENCE *** */
L70:
    *irc = -(i__ + *p);
    goto L200;

L80:
    i__ = -(*irc);
    if (i__ <= *p) {
	goto L100;
    }
    i__ -= *p;
    if (i__ > *p) {
	goto L90;
    }
    w[3] = *fx;
    h__ = w[5] * 2.;
    xih = w[6] + h__;
    *irc -= *p;
    goto L200;

/*    *** FINISH OFF-SIDE CENTRAL DIFFERENCE *** */

L90:
    i__ -= *p;
    g[i__] = (w[3] * 4. - *fx - w[4] * 3.) / w[5];
    *irc = i__;
    x[i__] = w[6];
    goto L20;

L100:
    h__ = -w[5];
    if (h__ > 0.) {
	goto L110;
    }
    w[3] = *fx;
    xih = w[6] + h__;
    goto L200;

L110:
    g[i__] = (w[3] - *fx) / (h__ * 2.);
    x[i__] = w[6];
    goto L20;

/*     ***  COMPUTE FORWARD DIFFERENCES IN VARIOUS CASES  *** */

L120:
    if (h__ >= axibar * .02) {
	h__ = h0 * axibar;
    }
    if (alphai * gi < 0.) {
	h__ = -h__;
    }
    goto L150;
L130:
    h__ = axibar;
    goto L150;
L140:
    h__ = h0 * axibar;

L150:
    hit = FALSE_;
L160:
    xih = xi + h__;
    if (h__ > 0.) {
	goto L170;
    }
    if (xih >= b[(i__ << 1) + 1]) {
	goto L200;
    }
    goto L180;
L170:
    if (xih <= b[(i__ << 1) + 2]) {
	goto L200;
    }
L180:
    if (hit) {
	goto L190;
    }
    hit = TRUE_;
    h__ = -h__;
    goto L160;

/*        *** ERROR RETURN... */
L190:
    *irc = i__ + *p;
    goto L230;

/*        *** RETURN FOR NEW FUNCTION VALUE... */
L200:
    x[i__] = xih;
    w[5] = h__;
    goto L999;

/*     ***  COMPUTE ACTUAL FORWARD DIFFERENCE  *** */

L210:
    g[*irc] = (*fx - w[4]) / w[5];
    x[*irc] = w[6];
    goto L20;

/*  ***  RESTORE FX AND INDICATE THAT G HAS BEEN COMPUTED  *** */

L220:
    *irc = 0;
L230:
    *fx = w[4];

L999:
    return 0;
/*  ***  LAST LINE OF DS3GRD FOLLOWS  *** */
} /* ds3grd_ */

#ifdef __cplusplus
	}
#endif
