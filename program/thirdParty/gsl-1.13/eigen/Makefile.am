noinst_LTLIBRARIES = libgsleigen.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_eigen.h
libgsleigen_la_SOURCES =  jacobi.c symm.c symmv.c nonsymm.c nonsymmv.c herm.c hermv.c gensymm.c gensymmv.c genherm.c genhermv.c gen.c genv.c sort.c francis.c schur.c

INCLUDES = -I$(top_srcdir)

noinst_HEADERS =  qrstep.c

TESTS = $(check_PROGRAMS)

test_LDADD = libgsleigen.la  ../test/libgsltest.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../utils/libutils.la ../rng/libgslrng.la ../sort/libgslsort.la

test_SOURCES = test.c


