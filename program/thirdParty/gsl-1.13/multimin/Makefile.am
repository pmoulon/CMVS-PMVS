noinst_LTLIBRARIES = libgslmultimin.la 

pkginclude_HEADERS = gsl_multimin.h

INCLUDES = -I$(top_srcdir)

libgslmultimin_la_SOURCES = fdfminimizer.c steepest_descent.c conjugate_fr.c conjugate_pr.c convergence.c diff.c vector_bfgs.c vector_bfgs2.c fminimizer.c simplex.c simplex2.c

noinst_HEADERS = directional_minimize.c linear_minimize.c linear_wrapper.c

check_PROGRAMS = test #demo

TESTS = $(check_PROGRAMS) 

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultimin.la ../min/libgslmin.la ../poly/libgslpoly.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#demo_SOURCES = demo.c 
#demo_LDADD = libgslmultimin.la ../min/libgslmin.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../linalg/libgsllinalg.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

