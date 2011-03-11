/* dgeqrf.f -- translated by f2c (version 20061008).
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
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* Subroutine */ int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;

    /* Local variables */
    integer i__, j, k, ib, nb, nt, nx, iws;
    extern doublereal sceil_(real *);
    integer nbmin, iinfo;
    extern /* Subroutine */ int dgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dlarfb_(char *, 
	     char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *), dlarft_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *);
    integer lbwork, llwork, lwkopt;
    logical lquery;


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     March 2008 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEQRF computes a QR factorization of a real M-by-N matrix A: */
/*  A = Q * R. */

/*  This is the left-looking Level 3 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, the elements on and above the diagonal of the array */
/*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/*          upper triangular if m >= n); the elements below the diagonal, */
/*          with the array TAU, represent the orthogonal matrix Q as a */
/*          product of min(m,n) elementary reflectors (see Further */
/*          Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N)) */
/*          The scalar factors of the elementary reflectors (see Further */
/*          Details). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */

/*          The dimension of the array WORK. The dimension can be divided into three parts. */

/*          1) The part for the triangular factor T. If the very last T is not bigger */
/*             than any of the rest, then this part is NB x ceiling(K/NB), otherwise, */
/*             NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T */

/*          2) The part for the very last T when T is bigger than any of the rest T. */
/*             The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB, */
/*             where K = min(M,N), NX is calculated by */
/*                   NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) ) */

/*          3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB) */

/*          So LWORK = part1 + part2 + part3 */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(k), where k = min(m,n). */

/*  Each H(i) has the form */

/*     H(i) = I - tau * v * v' */

/*  where tau is a real scalar, and v is a real vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/*  and tau in TAU(i). */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nbmin = 2;
    nx = 0;
    iws = *n;
    k = min(*m,*n);
    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1);
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1);
	nx = max(i__1,i__2);
    }

/*     Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.: */

/*            NB=3     2NB=6       K=10 */
/*            |        |           | */
/*      1--2--3--4--5--6--7--8--9--10 */
/*                  |     \________/ */
/*               K-NX=5      NT=4 */

/*     So here 4 x 4 is the last T stored in the workspace */

    r__1 = (real) (k - nx) / (real) nb;
    nt = k - sceil_(&r__1) * nb;

/*     optimal workspace = space for dlarfb + space for normal T's + space for the last T */

/* Computing MAX */
/* Computing MAX */
    i__3 = (*n - *m) * k, i__4 = (*n - *m) * nb;
/* Computing MAX */
    i__5 = k * nb, i__6 = nb * nb;
    i__1 = max(i__3,i__4), i__2 = max(i__5,i__6);
    llwork = max(i__1,i__2);
    r__1 = (real) llwork / (real) nb;
    llwork = sceil_(&r__1);
    if (nt > nb) {
	lbwork = k - nt;

/*         Optimal workspace for dlarfb = MAX(1,N)*NT */

	lwkopt = (lbwork + llwork) * nb;
	work[1] = (doublereal) (lwkopt + nt * nt);
    } else {
	r__1 = (real) k / (real) nb;
	lbwork = sceil_(&r__1) * nb;
	lwkopt = (lbwork + llwork - nb) * nb;
	work[1] = (doublereal) lwkopt;
    }

/*     Test the input arguments */

    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQRF", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (k == 0) {
	work[1] = 1.;
	return 0;
    }

    if (nb > 1 && nb < k) {
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    if (nt <= nb) {
		iws = (lbwork + llwork - nb) * nb;
	    } else {
		iws = (lbwork + llwork) * nb + nt * nt;
	    }
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		if (nt <= nb) {
		    nb = *lwork / (llwork + (lbwork - nb));
		} else {
		    nb = (*lwork - nt * nt) / (lbwork + llwork);
		}
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Update the current column using old T's */

	    i__3 = i__ - nb;
	    i__4 = nb;
	    for (j = 1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              Apply H' to A(J:M,I:I+IB-1) from the left */

		i__5 = *m - j + 1;
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__5, &
			ib, &nb, &a[j + j * a_dim1], lda, &work[j], &lbwork, &
			a[j + i__ * a_dim1], lda, &work[lbwork * nb + nt * nt 
			+ 1], &ib);
/* L20: */
	    }

/*           Compute the QR factorization of the current block */
/*           A(I:M,I:I+IB-1) */

	    i__4 = *m - i__ + 1;
	    dgeqr2_(&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    lbwork * nb + nt * nt + 1], &iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

		i__4 = *m - i__ + 1;
		dlarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[i__], &lbwork);

	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	if (i__ != 1) {
	    i__2 = i__ - nb;
	    i__1 = nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*                Apply H' to A(J:M,I:K) from the left */

		i__4 = *m - j + 1;
		i__3 = k - i__ + 1;
		i__5 = k - i__ + 1;
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__4, &
			i__3, &nb, &a[j + j * a_dim1], lda, &work[j], &lbwork, 
			 &a[j + i__ * a_dim1], lda, &work[lbwork * nb + nt * 
			nt + 1], &i__5);
/* L30: */
	    }
	    i__1 = *m - i__ + 1;
	    i__2 = k - i__ + 1;
	    dgeqr2_(&i__1, &i__2, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[lbwork * nb + nt * nt + 1], &iinfo);
	} else {

/*        Use unblocked code to factor the last or only block. */

	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__ + 1;
	    dgeqr2_(&i__1, &i__2, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
	}
    }

/*     Apply update to the column M+1:N when N > M */

    if (*m < *n && i__ != 1) {

/*         Form the last triangular factor of the block reflector */
/*         H = H(i) H(i+1) . . . H(i+ib-1) */

	if (nt <= nb) {
	    i__1 = *m - i__ + 1;
	    i__2 = k - i__ + 1;
	    dlarft_("Forward", "Columnwise", &i__1, &i__2, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], &work[i__], &lbwork);
	} else {
	    i__1 = *m - i__ + 1;
	    i__2 = k - i__ + 1;
	    dlarft_("Forward", "Columnwise", &i__1, &i__2, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], &work[lbwork * nb + 1], &nt);
	}

/*         Apply H' to A(1:M,M+1:N) from the left */

	i__1 = k - nx;
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__4 = k - j + 1;
	    ib = min(i__4,nb);
	    i__4 = *m - j + 1;
	    i__3 = *n - *m;
	    i__5 = *n - *m;
	    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__4, &
		    i__3, &ib, &a[j + j * a_dim1], lda, &work[j], &lbwork, &a[
		    j + (*m + 1) * a_dim1], lda, &work[lbwork * nb + nt * nt 
		    + 1], &i__5);
/* L40: */
	}
	if (nt <= nb) {
	    i__2 = *m - j + 1;
	    i__1 = *n - *m;
	    i__4 = k - j + 1;
	    i__3 = *n - *m;
	    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__2, &
		    i__1, &i__4, &a[j + j * a_dim1], lda, &work[j], &lbwork, &
		    a[j + (*m + 1) * a_dim1], lda, &work[lbwork * nb + nt * 
		    nt + 1], &i__3);
	} else {
	    i__2 = *m - j + 1;
	    i__1 = *n - *m;
	    i__4 = k - j + 1;
	    i__3 = *n - *m;
	    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__2, &
		    i__1, &i__4, &a[j + j * a_dim1], lda, &work[lbwork * nb + 
		    1], &nt, &a[j + (*m + 1) * a_dim1], lda, &work[lbwork * 
		    nb + nt * nt + 1], &i__3);
	}
    }
    work[1] = (doublereal) iws;
    return 0;

/*     End of DGEQRF */

} /* dgeqrf_ */
