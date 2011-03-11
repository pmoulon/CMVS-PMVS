noinst_LTLIBRARIES = libgslderiv.la

INCLUDES = -I$(top_srcdir)

libgslderiv_la_SOURCES = deriv.c

pkginclude_HEADERS = gsl_deriv.h

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test #demo

test_SOURCES = test.c
test_LDADD =  libgslderiv.la  ../vector/libgslvector.la ../block/libgslblock.la  ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

#demo_SOURCES = demo.c
#demo_LDADD =  libgslderiv.la  ../vector/libgslvector.la ../block/libgslblock.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la
