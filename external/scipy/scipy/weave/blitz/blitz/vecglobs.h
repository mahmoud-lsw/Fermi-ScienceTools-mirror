/***************************************************************************
 * blitz/vecglobs.h      Global vector functions
 *
 * $Id: vecglobs.h,v 1.1.2.2.8.1 2015/02/27 05:02:09 jasercio Exp $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

#ifndef BZ_VECGLOBS_H
#define BZ_VECGLOBS_H

#ifndef BZ_VECTOR_H
 #include <blitz/vector.h>
#endif

#ifndef BZ_NUMTRAIT_H
 #include <blitz/numtrait.h>
#endif

#ifndef BZ_PROMOTE_H
 #include <blitz/promote.h>
#endif

#ifndef BZ_EXTREMUM_H
 #include <blitz/extremum.h>
#endif


BZ_NAMESPACE(blitz)

BZ_NAMESPACE_END

#include <blitz/vecglobs.cc>

#endif // BZ_VECGLOBS_H
