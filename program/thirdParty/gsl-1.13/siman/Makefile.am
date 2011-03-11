## Process this file with automake to produce Makefile.in

check_PROGRAMS = test 
noinst_PROGRAMS = siman_tsp
noinst_LTLIBRARIES = libgslsiman.la

TESTS = $(check_PROGRAMS)
EXTRA_DIST = siman_test_driver.sh

test_SOURCES = test.c
test_LDADD = libgslsiman.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

siman_tsp_SOURCES = siman_tsp.c
siman_tsp_LDADD = libgslsiman.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../sys/libgslsys.la ../utils/libutils.la

CLEANFILES = siman_test.out

libgslsiman_la_SOURCES = siman.c
pkginclude_HEADERS = gsl_siman.h

INCLUDES = -I$(top_srcdir)
