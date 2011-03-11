noinst_LTLIBRARIES = libgslerr.la

pkginclude_HEADERS = gsl_errno.h gsl_message.h

libgslerr_la_SOURCES = error.c stream.c message.c strerror.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la
