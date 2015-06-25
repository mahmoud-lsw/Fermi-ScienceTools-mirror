# $Id: PSFMath.py,v 1.3 2008/09/05 19:40:49 elwinter Exp $

# Python module to provide numerical code for PyPSF.

# Most of these routines were adapted from:

# Numerical Recipes in C: The Art of Scientific Programming, 2nd ed.
# Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P.
# Cambridge University Press, 1992
# ISBN 0-521-43108-5

# In comments below, this citation is implied by the appearance of the
# string "(NRiC2e)".

# Reference values for some special functions are taken from:

# CRC Standard Mathematical Tables, 26th ed.
# Beyer, W.H., wd.
# CRC Press, 1981
# ISBN 0-8493-0626-4

# In comments below, this citation is implied by the appearance of the
# string "(SMT)".

# Reference values for some special functions are taken from:

# The Wolfram Functions Site
# http://functions.wolfram.com/
# Wolfram Research, Inc. 2008

# In comments below, this citation is implied by the appearance of the
# string "(WFS)".

#******************************************************************************

# Import external modules.

# Standard modules
from math import exp, floor, log

# Third-party modules

# Project modules

#******************************************************************************

# The chi-square probability function PChiSq(chi_sq, n_dof) computes
# the probability that the observed chi-square value for a correct
# model is less than the specified value chi_sq, using the specified
# number of degrees of freedom, n_dof. The algorithm is from
# NRiC.

# Example: If PChiSq(4.5, 10) = 0.95, then there is a 95% chance that
# the observed chi_sq value will be less than 4.5, assuming 10 degrees
# of freedom. So for any given pair (chi_sq0, n_dof), a high value of
# PChiSq means that there is a high probability that the observed
# chi_sq will be less than that value chi_sq0. So a high value of
# PChiSq for a low value of chi_sq0 indicates we have tightly
# constrained our model.

def PChiSq(chi_sq, n_dof):
    p = gammp(n_dof / 2.0, chi_sq / 2.0)
    return p

#------------------------------------------------------------------------------

# Compute the natural logarithm of the gamma function of xx.

# Based on gammln() (NRiC2e).

def gammln(xx):
    cof = [ 76.18009172947146,
           -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
             0.1208650973866179e-2,
            -0.5395239384953e-5 ]
    x = xx
    y = xx
    tmp = x + 5.5
    tmp -= (x + 0.5) * log(tmp)
    ser = 1.000000000190015
    for j in range(0, 6):
        y += 1
        ser += cof[j] / y
    return -tmp + log(2.5066282746310005 * ser / x)

#------------------------------------------------------------------------------

# Compute the floating-point factorial of n.

# Based on factrl() (NRiC2e).

def factrl(n):
    ntop = 4
    a = []
    for i in range(0, 33):
        a.append(0.0)
    a[0:5] = (1.0, 1.0, 2.0, 6.0, 24.0)
    assert(n >= 0)
    if n > 32:
        return exp(gammln(n + 1.0))
    while ntop < n:
        j = ntop
        ntop += 1
        a[ntop] = a[j] * ntop
    return a[n]

#------------------------------------------------------------------------------

# Compute the natural logarithm of the factorial of n.

# Based on factln() (NRiC2e).

def factln(n):
    a = []
    for i in range(0, 101):
        a.append(0.0)
    assert(n >= 0)
    if n <= 1:
        return 0.0
    if n <= 100:
        if a[n] > 0.0:
            return a[n]
        else:
            a[n] = gammln(n + 1.0)
            return a[n]
    else:
        return gammln(n + 1.0)

#------------------------------------------------------------------------------

# Compute the floating-point binomial coefficient of n and k.

# Based on bico() (NRiC2e).

def bico(n, k):
    return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)))

#------------------------------------------------------------------------------

# Compute the beta function of z and w.

# Based on beta() (NRiC2e).

def beta(z, w):
    return exp(gammln(z) + gammln(w) - gammln(z + w))

#------------------------------------------------------------------------------

# Compute the incomplete gamma function P(a,x) using a series
# representation. This algorithm converges rapidly for x <~ a + 1.

# Based on gser() (NRiC2e).

ITMAX = 100
EPS = 3.0e-7
FPMIN = 1.0e-30

def gser(gamser, a, x, gln):
    gln = gammln(a)
    if x <= 0.0:
        assert(x >= 0.0)
        gamser = 0.0
        return gamser
    else:
        ap = a
        delta = 1.0 / a
        sum = delta
        for n in range(1, ITMAX + 1):
            ap += 1
            delta *= x / float(ap)
            sum += delta
            if abs(delta) < abs(sum) * EPS:
                gamser = sum * exp(-x + a * log(x) - gln)
                return gamser
        assert(0)

#------------------------------------------------------------------------------

# Compute the incomplete gamma function Q(a,x) = 1 - P(a,x) using a
# continued fraction representation. This algorithm converges rapidly
# for x >~ a + 1.

# Based on gcf() (NRiC2e).

def gcf(gammcf, a, x, gln):
    gln = gammln(a)
    b = x + 1.0 - a
    c = 1.0 / FPMIN
    d = 1.0 / b
    h = d
    for i in range(1, ITMAX + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < FPMIN:
            d = FPMIN
        c = b + an / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    assert(i <= ITMAX)
    gammcf = exp(-x + a * log(x) - gln) * h
    return gammcf

#------------------------------------------------------------------------------

# Compute the incomplete gamma function P(a, x).

# Based on gammp() (NRiC2e).

def gammp(a, x):
    assert(x >= 0.0 and a > 0.0)
    gln = 0.0
    if x < a + 1.0:
        gamser = 0.0
        gamser = gser(gamser, a, x, gln)
        return gamser
    else:
        gammcf = 0.0
        gammcf = gcf(gammcf, a, x, gln)
        return 1.0 - gammcf

#------------------------------------------------------------------------------

# Compute the incomplete gamma function Q(a, x) = 1 - P(a, x).

# Based on gammq() (NRiC2e).

def gammq(a, x):
    assert(x >= 0.0 and a > 0.0)
    if x < a + 1.0:
        gamser = 0.0
        gamser = gser(gamser, a, x, gln)
        return 1.0 - gamser
    else:
        gammcf = 0.0
        gammcf = gcf(gammcf, a, x, gln)
        return gammcf

#******************************************************************************

# Self-test code.

# If this code generates any output, an error has been found.

# Absolute tolerance when comparing two floating-point numbers for
# equality.
ABSOLUTE_TOLERANCE = 1.0e-5

if __name__ == '__main__':

    from math import exp

    #--------------------------------------------------------------------------

    # Test the ln(gamma) function, noting that:
    # n! = gamma(n+1), where n is an integer
    for n in range(0, 10):
        fac = 1
        for i in range(n, 1, -1):
            fac *= i
        gamma_np1 = int(round(exp(gammln(n + 1))))
        assert(fac == gamma_np1)

    # Test the ln(gamma) function, noting that:
    # gamma(z+1) = z * gamma(z), where z is a real number.
    for i in range(0, 10):
        z = 0.5 + i
        gamma_z = exp(gammln(z))
        gamma_zp1 = exp(gammln(z + 1.0))
        assert(abs(gamma_zp1 - z * gamma_z) <= ABSOLUTE_TOLERANCE)

    # Test a few other values of the ln(gamma) function, taken from
    # standard mathematical tables (SMT).
    z_test = [ 1.01, 1.23, 1.45, 1.56, 1.78, 1.90 ]
    gamma_z_test = [ 0.99433, 0.91075, 0.88566, 0.88964, 0.92623, 0.96177 ]
    for i in range(0, len(z_test)):
        gamma_z = exp(gammln(z_test[i]))
        assert(abs(gamma_z - gamma_z_test[i]) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the factorial function.
    for n in range(0, 11):
        fac1 = 1
        for i in range(n, 1, -1):
            fac1 *= i
        fac2 = exp(gammln(n + 1))
        assert(abs(fac1 - fac2) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the ln(factorial) function.
    for n in range(0, 100):
        fac = 1
        for i in range(n, 1, -1):
            fac *= i
        ln_fac1 = log(fac)
        ln_fac2 = factln(n)
        assert(abs(ln_fac1 - ln_fac2) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the bico() function.
    bc1 = ( (),
            (1, 1),
            (1, 2, 1),
            (1, 3, 3, 1),
            (1, 4, 6, 4, 1) )
    for n in range(1, 5):
        for k in range(0, n + 1):
            bc2 = int(bico(n, k))
            assert(bc2 == bc1[n][k])

    #--------------------------------------------------------------------------

    # Test the beta function. Reference values taken from (WFS).
    beta1 = ((),
             (None, 1.00000, 0.50000, 0.33333),
             (None, 0.50000, 0.16667, 0.08333),
             (None, 0.33333, 0.08333, 0.03333))
    for z in range(1, 4):
        for w in range(1, 4):
            beta2 = beta(z, w)
            assert(abs(beta1[z][w] - beta2) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the computation of the incomplete gamma function P(a,x)
    # using the series representation. Values taken from (WFS).

    # The (WFS) computes only the numerator of Q(a,x), so list that
    # separately as Qwfs(a,x), then compute Q1(a,x) and P1(a,x).
    Qwfs = ((),
            (1.00000, 0.367879, 0.135335, 0.049787),
            (1.00000, 0.735759, 0.406006, 0.199148),
            (2.00000, 1.83940 , 1.35335 , 0.846380)
            )
    Q1 = [[],
          [None, None, None, None],
          [None, None, None, None],
          [None, None, None, None]
          ]
    P1 = [[],
          [None, None, None, None],
          [None, None, None, None],
          [None, None, None, None]
          ]
    for a in range(1, 4):
        gama = exp(gammln(a))
        for x in range(0, 4):
            Q1[a][x] = Qwfs[a][x] / gama
            P1[a][x] = 1.0 - Q1[a][x]
    gamser = 0.0
    gln = 0.0
    for a in range(1, 4):
        for x in range(0, a + 1):
            P2 = gser(gamser, a, x, gln)
            assert(abs(P2 - P1[a][x]) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the computation of the incomplete gamma function Q(a,x)
    # using the continued fraction representation.
    gammcf = 0.0
    for a in range(1, 4):
        for x in range(a + 1, 4):
            Q2 = gcf(gammcf, a, x, gln)
            assert(abs(Q2 - Q1[a][x]) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the computation of the incomplete gamma function P(a,x)
    # using the dynamically-determined representation.
    for a in range(1, 4):
        for x in range(0, 4):
            P2 = gammp(a, x)
            assert(abs(P2 - P1[a][x]) <= ABSOLUTE_TOLERANCE)

    #--------------------------------------------------------------------------

    # Test the computation of the incomplete gamma function Q(a,x)
    # using the dynamically-determined representation.
    for a in range(1, 4):
        for x in range(0, 4):
            Q2 = gammq(a, x)
            assert(abs(Q2 - Q1[a][x]) <= ABSOLUTE_TOLERANCE)
