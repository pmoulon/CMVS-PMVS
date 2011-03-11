noinst_LTLIBRARIES = libgsllinalg.la 

pkginclude_HEADERS = gsl_linalg.h

INCLUDES = -I$(top_srcdir)

libgsllinalg_la_SOURCES = multiply.c exponential.c tridiag.c tridiag.h lu.c luc.c hh.c qr.c qrpt.c lq.c ptlq.c svd.c householder.c householdercomplex.c hessenberg.c hesstri.c cholesky.c choleskyc.c symmtd.c hermtd.c bidiag.c balance.c balancemat.c

noinst_HEADERS =  givens.c apply_givens.c svdstep.c tridiag.h 

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_LDADD = libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c


