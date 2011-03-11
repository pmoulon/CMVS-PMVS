noinst_LTLIBRARIES = libgslrng.la 

pkginclude_HEADERS = gsl_rng.h

INCLUDES = -I$(top_srcdir)

libgslrng_la_SOURCES = borosh13.c cmrg.c coveyou.c default.c file.c fishman18.c fishman20.c fishman2x.c gfsr4.c knuthran2.c knuthran.c knuthran2002.c lecuyer21.c minstd.c mrg.c mt.c r250.c ran0.c ran1.c ran2.c ran3.c rand48.c rand.c random.c randu.c ranf.c ranlux.c ranlxd.c ranlxs.c ranmar.c rng.c slatec.c taus.c taus113.c transputer.c tt.c types.c uni32.c uni.c vax.c waterman14.c zuf.c inline.c

CLEANFILES = test.dat

noinst_HEADERS = schrage.c

test_SOURCES = test.c
test_LDADD = libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test

# benchmark_SOURCES = benchmark.c 
# benchmark_LDADD = libgslrng.la ../err/libgslerr.la ../utils/libutils.la
# rng_dump_SOURCES = rng-dump.c 
# rng_dump_LDADD = libgslrng.la ../err/libgslerr.la ../utils/libutils.la
