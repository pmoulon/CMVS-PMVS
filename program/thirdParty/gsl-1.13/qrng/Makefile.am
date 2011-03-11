noinst_LTLIBRARIES = libgslqrng.la 

pkginclude_HEADERS = gsl_qrng.h

INCLUDES = -I$(top_srcdir)

libgslqrng_la_SOURCES =	gsl_qrng.h qrng.c niederreiter-2.c sobol.c halton.c reversehalton.c inline.c

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslqrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

