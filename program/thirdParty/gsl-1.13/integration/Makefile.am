noinst_LTLIBRARIES = libgslintegration.la

INCLUDES = -I$(top_srcdir)

libgslintegration_la_SOURCES = qk15.c qk21.c qk31.c qk41.c qk51.c qk61.c qk.c qng.c qng.h qag.c	qags.c qagp.c workspace.c qcheb.c qawc.c qmomo.c qaws.c	qmomof.c qawo.c	qawf.c 

pkginclude_HEADERS = gsl_integration.h
noinst_HEADERS = qpsrt.c qpsrt2.c qelg.c qc25c.c qc25s.c qc25f.c ptsort.c util.c err.c positivity.c append.c initialise.c set_initial.c reset.c

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test

test_SOURCES = test.c tests.c tests.h
test_LDADD = libgslintegration.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la
