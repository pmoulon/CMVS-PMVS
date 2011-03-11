noinst_LTLIBRARIES = libgslmonte.la 

libgslmonte_la_SOURCES =  miser.c plain.c gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h gsl_monte.h vegas.c

pkginclude_HEADERS = gsl_monte.h gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test #demo

test_SOURCES = test.c 
test_LDADD = libgslmonte.la ../rng/libgslrng.la  ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

noinst_HEADERS = test_main.c

#demo_SOURCES= demo.c
#demo_LDADD = libgslmonte.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la


