noinst_LTLIBRARIES = libgslfft.la 

pkginclude_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_dft_complex.h gsl_dft_complex_float.h gsl_fft_complex_float.h gsl_fft_halfcomplex_float.h gsl_fft_real_float.h

INCLUDES = -I$(top_srcdir)

libgslfft_la_SOURCES =  dft.c fft.c

noinst_HEADERS = c_pass.h hc_pass.h real_pass.h signals.h signals_source.c c_main.c c_init.c c_pass_2.c c_pass_3.c c_pass_4.c c_pass_5.c c_pass_6.c c_pass_7.c c_pass_n.c c_radix2.c bitreverse.c bitreverse.h factorize.c factorize.h hc_init.c hc_pass_2.c hc_pass_3.c hc_pass_4.c hc_pass_5.c hc_pass_n.c hc_radix2.c hc_unpack.c real_init.c real_pass_2.c real_pass_3.c real_pass_4.c real_pass_5.c real_pass_n.c real_radix2.c real_unpack.c compare.h compare_source.c dft_source.c hc_main.c real_main.c test_complex_source.c test_real_source.c test_trap_source.c urand.c complex_internal.h

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c signals.c

test_LDADD = libgslfft.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#errs_LDADD = libgslfft.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la
#benchmark_LDADD = libgslfft.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

