/* zgetrf.f -- translated by f2c (version 20061008).
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

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {-1.,0.};
static integer c__1 = 1;

/* Subroutine */ int zgetrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, ipivstart, jpivstart, jp;
    doublecomplex tmp;
    integer kcols;
    doublereal sfmin;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *, 
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);
    integer nstep;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *, 
	     doublecomplex *, integer *);
    integer kahead;
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    doublereal pivmag;
    integer npived;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    integer kstart, ntopiv;
    extern /* Subroutine */ int zlaswp_(integer *, doublecomplex *, integer *, 
	     integer *, integer *, integer *, integer *);


/*  -- LAPACK routine (version 3.X) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     May 2008 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZGETRF computes an LU factorization of a general M-by-N matrix A */
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

/*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
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
	xerbla_("ZGETRF", &i__1);
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
	jp = j - 1 + izamax_(&i__2, &a[j + j * a_dim1], &c__1);
	ipiv[j] = jp;
/*        Permute just this column. */
	if (jp != j) {
	    i__2 = j + j * a_dim1;
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
	    i__2 = j + j * a_dim1;
	    i__3 = jp + j * a_dim1;
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	    i__2 = jp + j * a_dim1;
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
	}
/*        Apply pending permutations to L */
	ntopiv = 1;
	ipivstart = j;
	jpivstart = j - ntopiv;
	while(ntopiv < kahead) {
	    zlaswp_(&ntopiv, &a[jpivstart * a_dim1 + 1], lda, &ipivstart, &j, 
		    &ipiv[1], &c__1);
	    ipivstart -= ntopiv;
	    ntopiv <<= 1;
	    jpivstart -= ntopiv;
	}
/*        Permute U block to match L */
	zlaswp_(&kcols, &a[(j + 1) * a_dim1 + 1], lda, &kstart, &j, &ipiv[1], 
		&c__1);
/*        Factor the current column */
	pivmag = z_abs(&a[j + j * a_dim1]);
	if (pivmag != 0. && ! disnan_(&pivmag)) {
	    if (pivmag >= sfmin) {
		i__2 = *m - j;
		z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
		zscal_(&i__2, &z__1, &a[j + 1 + j * a_dim1], &c__1);
	    } else {
		i__2 = *m - j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = j + i__ + j * a_dim1;
		    z_div(&z__1, &a[j + i__ + j * a_dim1], &a[j + j * a_dim1])
			    ;
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
		}
	    }
	} else if (pivmag == 0. && *info == 0) {
	    *info = j;
	}
/*        Solve for U block. */
	ztrsm_("Left", "Lower", "No transpose", "Unit", &kahead, &kcols, &
		c_b1, &a[kstart + kstart * a_dim1], lda, &a[kstart + (j + 1) *
		 a_dim1], lda);
/*        Schur complement. */
	i__2 = *m - j;
	zgemm_("No transpose", "No transpose", &i__2, &kcols, &kahead, &c_b2, 
		&a[j + 1 + kstart * a_dim1], lda, &a[kstart + (j + 1) * 
		a_dim1], lda, &c_b1, &a[j + 1 + (j + 1) * a_dim1], lda);
    }
/*     Handle pivot permutations on the way out of the recursion */
    npived = nstep & -nstep;
    j = nstep - npived;
    while(j > 0) {
	ntopiv = j & -j;
	i__1 = j + 1;
	zlaswp_(&ntopiv, &a[(j - ntopiv + 1) * a_dim1 + 1], lda, &i__1, &
		nstep, &ipiv[1], &c__1);
	j -= ntopiv;
    }
/*     If short and wide, handle the rest of the columns. */
    if (*m < *n) {
	i__1 = *n - *m;
	zlaswp_(&i__1, &a[(*m + kcols + 1) * a_dim1 + 1], lda, &c__1, m, &
		ipiv[1], &c__1);
	i__1 = *n - *m;
	ztrsm_("Left", "Lower", "No transpose", "Unit", m, &i__1, &c_b1, &a[
		a_offset], lda, &a[(*m + kcols + 1) * a_dim1 + 1], lda);
    }
    return 0;

/*     End of ZGETRF */

} /* zgetrf_ */
