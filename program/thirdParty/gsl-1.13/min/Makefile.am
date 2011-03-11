noinst_LTLIBRARIES = libgslmin.la 

pkginclude_HEADERS = gsl_min.h

noinst_HEADERS = min.h

INCLUDES = -I$(top_srcdir)

libgslmin_la_SOURCES = fsolver.c golden.c brent.c convergence.c bracketing.c quad_golden.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c test_funcs.c test.h
test_LDADD = libgslmin.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la
