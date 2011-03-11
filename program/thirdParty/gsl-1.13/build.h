/* Compile subsequent inline functions as static functions */

#ifdef __GSL_BUILD_H__
#error build.h must not be included multiple times
#endif

#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif     

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else 
#error must be called with COMPILE_INLINE_STATIC
#endif
