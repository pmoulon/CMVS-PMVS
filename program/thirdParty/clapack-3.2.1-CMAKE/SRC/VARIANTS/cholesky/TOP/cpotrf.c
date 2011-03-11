/* cpotrf.f -- translated by f2c (version 20061008).
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
static real c_b18 = -1.f;
static real c_b19 = 1.f;

/* Subroutine */ int cpotrf_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer j, jb, nb;
    extern /* Subroutine */ int cherk_(char *, char *, integer *, integer *, 
	    real *, complex *, integer *, real *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, complex *, integer *, complex *, 
	    integer *);
    logical upper;
    extern /* Subroutine */ int cpotf2_(char *, integer *, complex *, integer 
	    *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *);


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     March 2008 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CPOTRF computes the Cholesky factorization of a real symmetric */
/*  positive definite matrix A. */

/*  The factorization has the form */
/*     A = U**H * U,  if UPLO = 'U', or */
/*     A = L  * L**H,  if UPLO = 'L', */
/*  where U is an upper triangular matrix and L is lower triangular. */

/*  This is the top-looking block version of the algorithm, calling Level 3 BLAS. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  Upper triangle of A is stored; */
/*          = 'L':  Lower triangle of A is stored. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA,N) */
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/*          N-by-N upper triangular part of A contains the upper */
/*          triangular part of the matrix A, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading N-by-N lower triangular part of A contains the lower */
/*          triangular part of the matrix A, and the strictly upper */
/*          triangular part of A is not referenced. */

/*          On exit, if INFO = 0, the factor U or L from the Cholesky */
/*          factorization A = U**H*U or A = L*L**H. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, the leading minor of order i is not */
/*                positive definite, and the factorization could not be */
/*                completed. */

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

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CPOTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "CPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

	cpotf2_(uplo, n, &a[a_offset], lda, info);
    } else {

/*        Use blocked code. */

	if (upper) {

/*           Compute the Cholesky factorization A = U'*U. */

	    i__1 = *n;
	    i__2 = nb;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = min(i__3,i__4);

/*              Compute the current block. */

		i__3 = j - 1;
		ctrsm_("Left", "Upper", "Conjugate Transpose", "Non-unit", &
			i__3, &jb, &c_b1, &a[a_dim1 + 1], lda, &a[j * a_dim1 
			+ 1], lda);
		i__3 = j - 1;
		cherk_("Upper", "Conjugate Transpose", &jb, &i__3, &c_b18, &a[
			j * a_dim1 + 1], lda, &c_b19, &a[j + j * a_dim1], lda);

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

		cpotf2_("Upper", &jb, &a[j + j * a_dim1], lda, info);
		if (*info != 0) {
		    goto L30;
		}
/* L10: */
	    }

	} else {

/*           Compute the Cholesky factorization A = L*L'. */

	    i__2 = *n;
	    i__1 = nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = min(i__3,i__4);

/*              Compute the current block. */

		i__3 = j - 1;
		ctrsm_("Right", "Lower", "Conjugate Transpose", "Non-unit", &
			jb, &i__3, &c_b1, &a[a_dim1 + 1], lda, &a[j + a_dim1], 
			 lda);
		i__3 = j - 1;
		cherk_("Lower", "No Transpose", &jb, &i__3, &c_b18, &a[j + 
			a_dim1], lda, &c_b19, &a[j + j * a_dim1], lda);

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

		cpotf2_("Lower", &jb, &a[j + j * a_dim1], lda, info);
		if (*info != 0) {
		    goto L30;
		}
/* L20: */
	    }
	}
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

/*     End of CPOTRF */

} /* cpotrf_ */
