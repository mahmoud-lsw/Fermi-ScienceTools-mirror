/*
 * $Id: vecpickio.cc,v 1.1.2.2.8.1 2015/02/27 05:02:09 jasercio Exp $
 *
 * Copyright (C) 1997 Todd Veldhuizen <tveldhui@oonumerics.org>
 * All rights reserved.  Please see <blitz/blitz.h> for terms and
 * conditions of use.
 *
 */

#ifndef BZ_VECPICKIO_CC
#define BZ_VECPICKIO_CC

#ifndef BZ_VECPICK_H
 #error <blitz/vecpickio.cc> must be included via <blitz/vecpick.h>
#endif // BZ_VECPICK_H

BZ_NAMESPACE(blitz)

template<typename P_numtype>
ostream& operator<<(ostream& os, const VectorPick<P_numtype>& x)
{
    Vector<P_numtype> y(x.length());
    y = x;
    os << y;
    return os;
}

BZ_NAMESPACE_END

#endif // BZ_VECPICKIO_CC
