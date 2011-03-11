noinst_LTLIBRARIES = libgslinterpolation.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_interp.h gsl_spline.h

libgslinterpolation_la_SOURCES = accel.c akima.c cspline.c interp.c linear.c integ_eval.h spline.c poly.c inline.c

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslinterpolation.la ../poly/libgslpoly.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../cblas/libgslcblas.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c

