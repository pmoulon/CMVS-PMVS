pkginclude_HEADERS = gsl_const.h gsl_const_cgs.h gsl_const_mks.h gsl_const_cgsm.h gsl_const_mksa.h gsl_const_num.h

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c

test_LDADD = ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la




