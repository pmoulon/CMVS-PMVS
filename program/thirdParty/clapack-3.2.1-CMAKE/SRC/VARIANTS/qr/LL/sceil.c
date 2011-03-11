/* sceil.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "blaswrap.h"

doublereal sceil_(real *a)
{
    /* System generated locals */
    real ret_val;


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     June 2008 */

/*     .. Scalar Arguments ..* */
/*     .. */

/*  ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements ..* */

    if (*a - (integer) (*a) == 0.f) {
	ret_val = *a;
    } else if (*a > 0.f) {
	ret_val = (real) ((integer) (*a) + 1);
    } else {
	ret_val = (real) ((integer) (*a));
    }
    return ret_val;

} /* sceil_ */
