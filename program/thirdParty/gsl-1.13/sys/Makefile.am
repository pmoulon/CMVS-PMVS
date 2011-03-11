noinst_LTLIBRARIES = libgslsys.la 

pkginclude_HEADERS = gsl_sys.h

libgslsys_la_SOURCES = minmax.c prec.c hypot.c log1p.c expm1.c coerce.c invhyp.c pow_int.c infnan.c fdiv.c fcmp.c ldfrexp.c

INCLUDES = -I$(top_srcdir)

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)
test_SOURCES = test.c
test_LDADD = libgslsys.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la libgslsys.la ../utils/libutils.la
