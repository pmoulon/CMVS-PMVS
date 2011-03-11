/* dgetrf.f -- translated by f2c (version 20061008).
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = 1.;
static doublereal c_b15 = -1.;

/* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    integer i__, j, ipivstart, jpivstart, jp;
    doublereal tmp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
, doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    integer kcols;
    doublereal sfmin;
    integer nstep;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer kahead;
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    integer npived;
    extern /* Subroutine */ int dlaswp_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *);
    integer kstart, ntopiv;


/*  -- LAPACK routine (version 3.X) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     May 2008 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETRF computes an LU factorization of a general M-by-N matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This code implements an iterative version of Sivan Toledo's recursive */
/*  LU algorithm[1].  For square matrices, this iterative versions should */
/*  be within a factor of two of the optimum number of memory transfers. */

/*  The pattern is as follows, with the large blocks of U being updated */
/*  in one call to DTRSM, and the dotted lines denoting sections that */
/*  have had all pending permutations applied: */

/*   1 2 3 4 5 6 7 8 */
/*  +-+-+---+-------+------ */
/*  | |1|   |       | */
/*  |.+-+ 2 |       | */
/*  | | |   |       | */
/*  |.|.+-+-+   4   | */
/*  | | | |1|       | */
/*  | | |.+-+       | */
/*  | | | | |       | */
/*  |.|.|.|.+-+-+---+  8 */
/*  | | | | | |1|   | */
/*  | | | | |.+-+ 2 | */
/*  | | | | | | |   | */
/*  | | | | |.|.+-+-+ */
/*  | | | | | | | |1| */
/*  | | | | | | |.+-+ */
/*  | | | | | | | | | */
/*  |.|.|.|.|.|.|.|.+----- */
/*  | | | | | | | | | */

/*  The 1-2-1-4-1-2-1-8-... pattern is the position of the last 1 bit in */
/*  the binary expansion of the current column.  Each Schur update is */
/*  applied as soon as the necessary portion of U is available. */

/*  [1] Toledo, S. 1997. Locality of Reference in LU Decomposition with */
/*  Partial Pivoting. SIAM J. Matrix Anal. Appl. 18, 4 (Oct. 1997), */
/*  1065-1081. http://dx.doi.org/10.1137/S0895479896297744 */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix to be factored. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
/*                has been completed, but the factor U is exactly */
/*                singular, and division by zero will occur if it is used */
/*                to solve a system of equations. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Compute machine safe minimum */

    sfmin = dlamch_("S");

    nstep = min(*m,*n);
    i__1 = nstep;
    for (j = 1; j <= i__1; ++j) {
	kahead = j & -j;
	kstart = j + 1 - kahead;
/* Computing MIN */
	i__2 = kahead, i__3 = *m - j;
	kcols = min(i__2,i__3);

/*        Find pivot. */

	i__2 = *m - j + 1;
	jp = j - 1 + idamax_(&i__2, &a[j + j * a_dim1], &c__1);
	ipiv[j] = jp;
/*        Permute just this column. */
	if (jp != j) {
	    tmp = a[j + j * a_dim1];
	    a[j + j * a_dim1] = a[jp + j * a_dim1];
	    a[jp + j * a_dim1] = tmp;
	}
/*        Apply pending permutations to L */
	ntopiv = 1;
	ipivstart = j;
	jpivstart = j - ntopiv;
	while(ntopiv < kahead) {
	    dlaswp_(&ntopiv, &a[jpivstart * a_dim1 + 1], lda, &ipivstart, &j, 
		    &ipiv[1], &c__1);
	    ipivstart -= ntopiv;
	    ntopiv <<= 1;
	    jpivstart -= ntopiv;
	}
/*        Permute U block to match L */
	dlaswp_(&kcols, &a[(j + 1) * a_dim1 + 1], lda, &kstart, &j, &ipiv[1], 
		&c__1);
/*        Factor the current column */
	if (a[j + j * a_dim1] != 0. && ! disnan_(&a[j + j * a_dim1])) {
	    if ((d__1 = a[j + j * a_dim1], abs(d__1)) >= sfmin) {
		i__2 = *m - j;
		d__1 = 1. / a[j + j * a_dim1];
		dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    } else {
		i__2 = *m - j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
		}
	    }
	} else if (a[j + j * a_dim1] == 0. && *info == 0) {
	    *info = j;
	}
/*        Solve for U block. */
	dtrsm_("Left", "Lower", "No transpose", "Unit", &kahead, &kcols, &
		c_b12, &a[kstart + kstart * a_dim1], lda, &a[kstart + (j + 1) 
		* a_dim1], lda);
/*        Schur complement. */
	i__2 = *m - j;
	dgemm_("No transpose", "No transpose", &i__2, &kcols, &kahead, &c_b15, 
		 &a[j + 1 + kstart * a_dim1], lda, &a[kstart + (j + 1) * 
		a_dim1], lda, &c_b12, &a[j + 1 + (j + 1) * a_dim1], lda);
    }
/*     Handle pivot permutations on the way out of the recursion */
    npived = nstep & -nstep;
    j = nstep - npived;
    while(j > 0) {
	ntopiv = j & -j;
	i__1 = j + 1;
	dlaswp_(&ntopiv, &a[(j - ntopiv + 1) * a_dim1 + 1], lda, &i__1, &
		nstep, &ipiv[1], &c__1);
	j -= ntopiv;
    }
/*     If short and wide, handle the rest of the columns. */
    if (*m < *n) {
	i__1 = *n - *m;
	dlaswp_(&i__1, &a[(*m + kcols + 1) * a_dim1 + 1], lda, &c__1, m, &
		ipiv[1], &c__1);
	i__1 = *n - *m;
	dtrsm_("Left", "Lower", "No transpose", "Unit", m, &i__1, &c_b12, &a[
		a_offset], lda, &a[(*m + kcols + 1) * a_dim1 + 1], lda);
    }
    return 0;

/*     End of DGETRF */

} /* dgetrf_ */
