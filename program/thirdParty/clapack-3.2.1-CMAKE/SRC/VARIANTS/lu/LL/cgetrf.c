/* cgetrf.f -- translated by f2c (version 20061008).
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

static complex c_b1 = {1.f,0.f};
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int cgetrf_(integer *m, integer *n, complex *a, integer *lda, 
	 integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1;

    /* Local variables */
    integer i__, j, k, jb, nb;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, integer *);
    integer iinfo;
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, complex *, integer *, complex *, 
	    integer *), cgetf2_(integer *, 
	    integer *, complex *, integer *, integer *, integer *), xerbla_(
	    char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *);
    extern /* Subroutine */ int claswp_(integer *, complex *, integer *, 
	    integer *, integer *, integer *, integer *);


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     March 2008 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CGETRF computes an LU factorization of a general M-by-N matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This is the left-looking Level 3 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA,N) */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
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
	xerbla_("CGETRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "CGETRF", " ", m, n, &c_n1, &c_n1);
    if (nb <= 1 || nb >= min(*m,*n)) {

/*        Use unblocked code. */

	cgetf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

	i__1 = min(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = min(*m,*n) - j + 1;
	    jb = min(i__3,nb);


/*           Update before factoring the current panel */

	    i__3 = j - nb;
	    i__4 = nb;
	    for (k = 1; i__4 < 0 ? k >= i__3 : k <= i__3; k += i__4) {

/*              Apply interchanges to rows K:K+NB-1. */

		i__5 = k + nb - 1;
		claswp_(&jb, &a[j * a_dim1 + 1], lda, &k, &i__5, &ipiv[1], &
			c__1);

/*              Compute block row of U. */

		ctrsm_("Left", "Lower", "No transpose", "Unit", &nb, &jb, &
			c_b1, &a[k + k * a_dim1], lda, &a[k + j * a_dim1], 
			lda);

/*              Update trailing submatrix. */

		i__5 = *m - k - nb + 1;
		q__1.r = -1.f, q__1.i = -0.f;
		cgemm_("No transpose", "No transpose", &i__5, &jb, &nb, &q__1, 
			 &a[k + nb + k * a_dim1], lda, &a[k + j * a_dim1], 
			lda, &c_b1, &a[k + nb + j * a_dim1], lda);
/* L30: */
	    }

/*           Factor diagonal and subdiagonal blocks and test for exact */
/*           singularity. */

	    i__4 = *m - j + 1;
	    cgetf2_(&i__4, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__3 = *m, i__5 = j + jb - 1;
	    i__4 = min(i__3,i__5);
	    for (i__ = j; i__ <= i__4; ++i__) {
		ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
	    }

/* L20: */
	}

/*        Apply interchanges to the left-overs */

	i__2 = min(*m,*n);
	i__1 = nb;
	for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
	    i__4 = k - 1;
/* Computing MIN */
	    i__5 = k + nb - 1, i__6 = min(*m,*n);
	    i__3 = min(i__5,i__6);
	    claswp_(&i__4, &a[a_dim1 + 1], lda, &k, &i__3, &ipiv[1], &c__1);
/* L40: */
	}

/*        Apply update to the M+1:N columns when N > M */

	if (*n > *m) {
	    i__1 = *n - *m;
	    claswp_(&i__1, &a[(*m + 1) * a_dim1 + 1], lda, &c__1, m, &ipiv[1], 
		     &c__1);
	    i__1 = *m;
	    i__2 = nb;
	    for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
		i__4 = *m - k + 1;
		jb = min(i__4,nb);

		i__4 = *n - *m;
		ctrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__4, &
			c_b1, &a[k + k * a_dim1], lda, &a[k + (*m + 1) * 
			a_dim1], lda);

		if (k + nb <= *m) {
		    i__4 = *m - k - nb + 1;
		    i__3 = *n - *m;
		    q__1.r = -1.f, q__1.i = -0.f;
		    cgemm_("No transpose", "No transpose", &i__4, &i__3, &nb, 
			    &q__1, &a[k + nb + k * a_dim1], lda, &a[k + (*m + 
			    1) * a_dim1], lda, &c_b1, &a[k + nb + (*m + 1) * 
			    a_dim1], lda);
		}
/* L50: */
	    }
	}

    }
    return 0;

/*     End of CGETRF */

} /* cgetrf_ */
