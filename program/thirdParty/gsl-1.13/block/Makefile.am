noinst_LTLIBRARIES = libgslblock.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_block.h gsl_block_char.h gsl_block_complex_double.h gsl_block_complex_float.h gsl_block_complex_long_double.h gsl_block_double.h gsl_block_float.h gsl_block_int.h gsl_block_long.h gsl_block_long_double.h gsl_block_short.h gsl_block_uchar.h gsl_block_uint.h gsl_block_ulong.h gsl_block_ushort.h gsl_check_range.h

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = block_source.c init_source.c fprintf_source.c fwrite_source.c test_complex_source.c test_source.c test_io.c test_complex_io.c

libgslblock_la_SOURCES = init.c file.c block.c
