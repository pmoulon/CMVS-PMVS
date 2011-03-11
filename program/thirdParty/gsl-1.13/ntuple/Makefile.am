noinst_LTLIBRARIES = libgslntuple.la

pkginclude_HEADERS = gsl_ntuple.h

INCLUDES = -I$(top_srcdir)

libgslntuple_la_SOURCES = ntuple.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test #demo demo1

test_SOURCES = test.c
test_LDADD = libgslntuple.la ../histogram/libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#demo_SOURCES = demo.c
#demo_LDADD = libgslntuple.la ../histogram/libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#demo1_SOURCES = demo1.c
#demo1_LDADD = libgslntuple.la ../histogram/libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

CLEANFILES = test.dat
