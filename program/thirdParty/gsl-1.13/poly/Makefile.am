noinst_LTLIBRARIES = libgslpoly.la 

pkginclude_HEADERS = gsl_poly.h

INCLUDES = -I$(top_srcdir)

libgslpoly_la_SOURCES = dd.c eval.c solve_quadratic.c solve_cubic.c zsolve_quadratic.c zsolve_cubic.c zsolve.c zsolve_init.c deriv.c

noinst_HEADERS = balance.c companion.c qr.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslpoly.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

