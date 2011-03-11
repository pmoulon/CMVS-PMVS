noinst_LTLIBRARIES = libgslcombination.la 

pkginclude_HEADERS = gsl_combination.h

INCLUDES = -I$(top_srcdir)

libgslcombination_la_SOURCES = init.c file.c combination.c inline.c

noinst_HEADERS = 

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c

test_LDADD = libgslcombination.la ../vector/libgslvector.la  ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#noinst_PROGRAMS = demo
#demo_SOURCES = demo.c
#demo_LDADD = libgslcombination.la ../vector/libgslvector.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#CLEANFILES = test.txt test.dat


