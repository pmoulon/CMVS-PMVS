dnl Process this file with autoconf to produce a configure script.

AC_INIT([gsl],[1.13])
AC_CONFIG_SRCDIR(gsl_math.h)

AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

dnl Library versioning (C:R:A == current:revision:age)
dnl See the libtool manual for an explanation of the numbers
dnl
dnl gsl-1.0    libgsl 0:0:0  libgslcblas 0:0:0
dnl gsl-1.1    libgsl 1:0:1  libgslcblas 0:0:0
dnl gsl-1.1.1  libgsl 2:0:2  libgslcblas 0:0:0
dnl gsl-1.2    libgsl 3:0:3  libgslcblas 0:0:0
dnl gsl-1.3    libgsl 4:0:4  libgslcblas 0:0:0
dnl gsl-1.4    libgsl 5:0:5  libgslcblas 0:0:0
dnl gsl-1.5    libgsl 6:0:6  libgslcblas 0:0:0
dnl gsl-1.6    libgsl 7:0:7  libgslcblas 0:0:0
dnl gsl-1.7    libgsl 8:0:8  libgslcblas 0:0:0
dnl gsl-1.8    libgsl 9:0:9  libgslcblas 0:0:0
dnl gsl-1.9    libgsl 10:0:10 libgslcblas 0:0:0 
dnl gsl-1.10   libgsl 10:0:10 (*) libgslcblas 0:0:0 
dnl gsl-1.11   libgsl 12:0:12  libgslcblas 0:0:0 
dnl gsl-1.12   libgsl 13:0:13  libgslcblas 0:0:0 
dnl gsl-1.13   libgsl 14:0:14  libgslcblas 0:0:0 
dnl 
dnl (*) There was an error on this release.  Firstly, the versioning
dnl numbers were not updated.  Secondly, 2 functions were removed, but
dnl the age not reset--this should have been 11:0:0.  However these
dnl functions were not documented and are regarded as internal, so we
dnl will assume 11:0:11.
dnl
dnl How to update library version number
dnl ====================================
dnl 
dnl C: increment if the interface has additions, changes, removals.
dnl
dnl R: increment any time the source changes; set to 0 if you
dnl incremented CURRENT
dnl
dnl A: increment if any interfaces have been added; set to 0 if any
dnl interfaces have been removed. removal has precedence over adding,
dnl so set to 0 if both happened.

dnl
GSL_CURRENT=14
GSL_REVISION=0
GSL_AGE=14
dnl
CBLAS_CURRENT=0
CBLAS_REVISION=0
CBLAS_AGE=0

GSL_LT_VERSION="${GSL_CURRENT}:${GSL_REVISION}:${GSL_AGE}"
AC_SUBST(GSL_LT_VERSION)

GSL_LT_CBLAS_VERSION="${CBLAS_CURRENT}:${CBLAS_REVISION}:${CBLAS_AGE}"
AC_SUBST(GSL_LT_CBLAS_VERSION)

case "$VERSION" in
    *+)
        ;;
    *)
        AC_DEFINE(RELEASED,[],[Defined if this is an official release])
        ;;
esac

dnl things required by automake
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_LANG(C)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL

dnl Check compiler features
AC_TYPE_SIZE_T
dnl AC_C_CONST
AC_C_VOLATILE
AC_C_INLINE
AC_C_CHAR_UNSIGNED

GSL_CFLAGS="-I$includedir"
GSL_LIBS="-L$libdir -lgsl"

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

if test "$ac_cv_c_inline" != no ; then 
dnl Check for "extern inline", using a modified version of the test
dnl for AC_C_INLINE from acspecific.mt
dnl
   AC_CACHE_CHECK([for GNU-style extern inline], ac_cv_c_extern_inline,
   [ac_cv_c_extern_inline=no
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern $ac_cv_c_inline double foo(double x);
   extern $ac_cv_c_inline double foo(double x) { return x + 1.0 ; } ;
   double foo (double x) { return x + 1.0 ; };]], [[  foo(1.0)  ]])],[ac_cv_c_extern_inline="yes"],[])
   ])

   if test "$ac_cv_c_extern_inline" != no ; then
      AC_DEFINE(HAVE_INLINE,[1],[Define if you have inline])
   else
      AC_CACHE_CHECK([for C99-style inline], ac_cv_c_c99inline,
      [ac_cv_c_c99inline=no
      dnl next line is a necessary condition
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern inline void* foo() { foo(); return &foo ; };]], 
      [[  return foo() != 0 ]])],[ac_cv_c_c99inline="yes"],[])
      dnl but not sufficient, extern must work but inline on its own should not
      if test "$ac_cv_c_c99inline" != no ; then
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[inline void* foo() { foo(); return &foo ; };]], 
            [[  return foo() != 0 ]])],[],ac_cv_c_c99inline="no")
      fi
      ])
      if test "$ac_cv_c_c99inline" != no ; then
         AC_DEFINE(HAVE_INLINE,[1],[Define if you have inline])
         AC_DEFINE(HAVE_C99_INLINE,[1],[Define if you have inline with C99 behavior])
      fi
   fi
fi

dnl Checks for header files.
AC_CHECK_HEADERS(ieeefp.h)

dnl Checks for typedefs, structures, and compiler characteristics.

case "$host_os" in
     *mingw*)
	MINGW=true
	;;
esac

AM_CONDITIONAL(MINGW32_HOST, test "$MINGW" = "true")

dnl Checks for library functions.

dnl AC_FUNC_ALLOCA
AC_FUNC_VPRINTF

dnl strcasecmp, strerror, xmalloc, xrealloc, probably others should be added.
dnl removed strerror from this list, it's hardcoded in the err/ directory
dnl Any functions which appear in this list of functions should be provided
dnl in the utils/ directory
dnl xmalloc is not used, removed (bjg)
AC_REPLACE_FUNCS(memcpy memmove strdup strtol strtoul)

AC_CACHE_CHECK(for EXIT_SUCCESS and EXIT_FAILURE,
ac_cv_decl_exit_success_and_failure,
AC_EGREP_CPP(yes,
[
#include <stdlib.h>
#ifdef EXIT_SUCCESS
yes
#endif
], 
ac_cv_decl_exit_success_and_failure=yes,
ac_cv_decl_exit_success_and_failure=no)
)

if test "$ac_cv_decl_exit_success_and_failure" = yes ; then
  AC_DEFINE(HAVE_EXIT_SUCCESS_AND_FAILURE,1,[Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h])
fi ;

dnl Use alternate libm if specified by user

if test "x$LIBS" = "x" ; then
  AC_CHECK_LIB(m, cos)
fi

dnl Remember to put a definition in acconfig.h for each of these
AC_CHECK_DECLS(feenableexcept,,,[#define _GNU_SOURCE 1
#include <fenv.h>]) 
AC_CHECK_DECLS(fesettrapenable,,,[#define _GNU_SOURCE 1
#include <fenv.h>]) 
AC_CHECK_DECLS(hypot,,,[#include <math.h>]) 
AC_CHECK_DECLS(expm1,,,[#include <math.h>])
AC_CHECK_DECLS(acosh,,,[#include <math.h>])
AC_CHECK_DECLS(asinh,,,[#include <math.h>])
AC_CHECK_DECLS(atanh,,,[#include <math.h>])
AC_CHECK_DECLS(ldexp,,,[#include <math.h>])
AC_CHECK_DECLS(frexp,,,[#include <math.h>])
AC_CHECK_DECLS(isinf,,,[#include <math.h>])
AC_CHECK_DECLS(isfinite,,,[#include <math.h>])
AC_CHECK_DECLS(finite,,,[#include <math.h>
#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif])
AC_CHECK_DECLS(isnan,,,[#include <math.h>])

dnl OpenBSD has a broken implementation of log1p.
case "$host" in
    *-*-*openbsd*)
       AC_MSG_RESULT([avoiding OpenBSD system log1p - using gsl version])
       ;;
    *)
        AC_CHECK_DECLS(log1p,,,[#include <math.h>])
       ;;
esac

AC_CACHE_CHECK([for long double stdio], ac_cv_func_printf_longdouble,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include <stdio.h>
int main (void) 
{ 
const char * s = "5678.25"; long double x = 1.234 ; 
fprintf(stderr,"%Lg\n",x) ; 
sscanf(s, "%Lg", &x);
if (x == 5678.25) {exit (0);} else {exit(1); };
}]])],[ac_cv_func_printf_longdouble="yes"],[ac_cv_func_printf_longdouble="no"],[ac_cv_func_printf_longdouble="no"])])

if test "$ac_cv_func_printf_longdouble" != no; then
  AC_DEFINE(HAVE_PRINTF_LONGDOUBLE,1,[Define this if printf can handle %Lf for long double])
fi

AC_CACHE_CHECK([for extended floating point registers],ac_cv_c_extended_fp,
[case "$host" in
    *sparc*-*-*)
        ac_cv_c_extended_fp=no
        ;;     
    *powerpc*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *hppa*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *alpha*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *68k*-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    *86-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    x86_64-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    *) 
        ac_cv_c_extended_fp=unknown
        ;;
esac
])

if test $ac_cv_c_extended_fp != "no" ; then
    AC_DEFINE(HAVE_EXTENDED_PRECISION_REGISTERS,1,[Defined on architectures with excess floating-point precision])
fi

AC_CACHE_CHECK([for IEEE arithmetic interface type], ac_cv_c_ieee_interface,
[case "$host" in
    sparc-*-linux*) 
        ac_cv_c_ieee_interface=gnusparc
        ;;
    m68k-*-linux*) 
        ac_cv_c_ieee_interface=gnum68k
        ;;
    powerpc-*-linux*) 
        ac_cv_c_ieee_interface=gnuppc
        ;;
    *86-*-gnu | *86_64-*-gnu | *86-*-linux* | *86_64-*-linux*) 
        ac_cv_c_ieee_interface=gnux86
        ;;
    *-*-sunos4*) 
        ac_cv_c_ieee_interface=sunos4
        ;;
    *-*-solaris*) 
        ac_cv_c_ieee_interface=solaris
        ;;
    *-*-hpux11*) 
        ac_cv_c_ieee_interface=hpux11
        ;;
    *-*-hpux*) 
        ac_cv_c_ieee_interface=hpux
        ;;
    *-*-osf*) 
        ac_cv_c_ieee_interface=tru64
        ;;
    *-*-aix*) 
        ac_cv_c_ieee_interface=aix
        ;;
    *-*-irix*) 
        ac_cv_c_ieee_interface=irix
        ;;
    powerpc-*-*darwin*) 
        ac_cv_c_ieee_interface=darwin
        ;;
    *86-*-*darwin*) 
        ac_cv_c_ieee_interface=darwin86
        ;;
    *-*-*netbsd*) 
        ac_cv_c_ieee_interface=netbsd
        ;;
    *-*-*openbsd*)  
        ac_cv_c_ieee_interface=openbsd
        ;;
    *-*-*bsd*) 
        ac_cv_c_ieee_interface=freebsd
        ;;
    *-*-os2*)
        ac_cv_c_ieee_interface=os2emx
        ;;
    *)
        ac_cv_c_ieee_interface=unknown
        ;;
esac
])

if test "$ac_cv_c_ieee_interface" = "gnux86" ; then
    AC_CACHE_CHECK([for FPU_SETCW], ac_cv_c_fpu_setcw,
    [ac_cv_c_fpu_setcw=no
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <fpu_control.h>
#ifndef _FPU_SETCW
#include <i386/fpu_control.h>
#define _FPU_SETCW(cw) __setfpucw(cw)
#endif
]], [[ unsigned short mode = 0 ; _FPU_SETCW(mode); ]])],[ac_cv_c_fpu_setcw="yes"],[ac_cv_c_ieee_interface=unknown])
    ])
fi

if test "$ac_cv_c_ieee_interface" = "gnux86" ; then
    AC_CACHE_CHECK([for SSE extensions], ac_cv_c_fpu_sse,
    [ac_cv_c_fpu_sse=no
    AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <stdlib.h>
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (*&cw_sse))
]], [[ unsigned int mode = 0x1f80 ; _FPU_SETMXCSR(mode); exit(0); ]])],[ac_cv_c_fpu_sse="yes"],[ac_cv_c_fpu_sse="no"],[
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stdlib.h>
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (*&cw_sse))
]], [[ unsigned int mode = 0x1f80 ; _FPU_SETMXCSR(mode); exit(0); ]])],[ac_cv_c_fpu_sse="yes"],[ac_cv_c_fpu_sse="no"])
])])

    if test $ac_cv_c_fpu_sse = yes; then
        AC_DEFINE([HAVE_FPU_X86_SSE], 1,
                  [Define if x86 processor has sse extensions.])
   fi
fi

ac_tr_ieee_interface=HAVE_`echo $ac_cv_c_ieee_interface | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`_IEEE_INTERFACE
AC_DEFINE_UNQUOTED($ac_tr_ieee_interface,1,[IEEE Interface Type])

AC_SUBST(HAVE_GNUSPARC_IEEE_INTERFACE)
AC_SUBST(HAVE_GNUM68K_IEEE_INTERFACE)
AC_SUBST(HAVE_GNUPPC_IEEE_INTERFACE)
AC_SUBST(HAVE_GNUX86_IEEE_INTERFACE)
AC_SUBST(HAVE_SUNOS4_IEEE_INTERFACE)
AC_SUBST(HAVE_SOLARIS_IEEE_INTERFACE)
AC_SUBST(HAVE_HPUX11_IEEE_INTERFACE)
AC_SUBST(HAVE_HPUX_IEEE_INTERFACE)
AC_SUBST(HAVE_TRU64_IEEE_INTERFACE)
AC_SUBST(HAVE_IRIX_IEEE_INTERFACE)
AC_SUBST(HAVE_AIX_IEEE_INTERFACE)
AC_SUBST(HAVE_FREEBSD_IEEE_INTERFACE)
AC_SUBST(HAVE_OS2EMX_IEEE_INTERFACE)
AC_SUBST(HAVE_NETBSD_IEEE_INTERFACE)
AC_SUBST(HAVE_OPENBSD_IEEE_INTERFACE)
AC_SUBST(HAVE_DARWIN_IEEE_INTERFACE)
AC_SUBST(HAVE_DARWIN86_IEEE_INTERFACE)

dnl Check for IEEE control flags

save_cflags="$CFLAGS"
AC_CACHE_CHECK([for IEEE compiler flags], ac_cv_c_ieee_flags,
[
case "$host" in
    alpha*-*-*)
        if test X"$GCC" = Xyes ; then
            ieee_flags='-mieee -mfp-rounding-mode=d'
        else
            # This assumes Compaq's C compiler.
            ieee_flags='-ieee -fprm d'
        fi
        ;;
esac
if test X"$ieee_flags" != X ; then
  CFLAGS="$ieee_flags $CFLAGS"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[int foo;]])],[ac_cv_c_ieee_flags="$ieee_flags"],[ac_cv_c_ieee_flags="none"])
else
  ac_cv_c_ieee_flags="none"
fi])

if test "$ac_cv_c_ieee_flags" != "none" ; then
   CFLAGS="$ac_cv_c_ieee_flags $save_cflags"
else
   CFLAGS="$save_cflags"
fi

dnl Check IEEE comparisons, whether "x != x" is true for NaNs
dnl
AC_CACHE_CHECK([for IEEE comparisons], ac_cv_c_ieee_comparisons,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <math.h>
int main (void) 
{ 
   int status; double inf, nan;
   inf = exp(1.0e10);
   nan = inf / inf ;
   status = (nan == nan);
   exit (status);
}]])],[ac_cv_c_ieee_comparisons="yes"],[ac_cv_c_ieee_comparisons="no"],[ac_cv_c_ieee_comparisons="yes"])
])

if test "$ac_cv_c_ieee_comparisons" != no ; then
  AC_DEFINE(HAVE_IEEE_COMPARISONS,1,[Define this if IEEE comparisons work correctly (e.g. NaN != NaN)])
fi

dnl Check for IEEE denormalized arithmetic
dnl
AC_CACHE_CHECK([for IEEE denormalized values], ac_cv_c_ieee_denormals,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <math.h> 
int main (void) 
{ 
   int i, status; 
   volatile double z = 1e-308;
   for (i = 0; i < 5; i++) { z = z / 10.0 ; };        
   for (i = 0; i < 5; i++) { z = z * 10.0 ; };
   status = (z == 0.0);
   exit (status);
}]])],[ac_cv_c_ieee_denormals="yes"],[ac_cv_c_ieee_denormals="no"],[ac_cv_c_ieee_denormals="yes"])
])

if test "$ac_cv_c_ieee_denormals" != no ; then
  AC_DEFINE(HAVE_IEEE_DENORMALS,1,[Define this if IEEE denormalized numbers are available])
fi

AH_TEMPLATE([HIDE_INLINE_STATIC],[Define if you need to hide the static definitions of inline functions])

AH_BOTTOM([/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#if !HAVE_EXIT_SUCCESS_AND_FAILURE
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif])

AH_BOTTOM([/* Define one of these if you have a known IEEE arithmetic interface */
#undef HAVE_GNUSPARC_IEEE_INTERFACE
#undef HAVE_GNUM68K_IEEE_INTERFACE
#undef HAVE_GNUPPC_IEEE_INTERFACE
#undef HAVE_GNUX86_IEEE_INTERFACE
#undef HAVE_SUNOS4_IEEE_INTERFACE
#undef HAVE_SOLARIS_IEEE_INTERFACE
#undef HAVE_HPUX11_IEEE_INTERFACE
#undef HAVE_HPUX_IEEE_INTERFACE
#undef HAVE_TRU64_IEEE_INTERFACE
#undef HAVE_IRIX_IEEE_INTERFACE
#undef HAVE_AIX_IEEE_INTERFACE
#undef HAVE_FREEBSD_IEEE_INTERFACE
#undef HAVE_OS2EMX_IEEE_INTERFACE
#undef HAVE_NETBSD_IEEE_INTERFACE
#undef HAVE_OPENBSD_IEEE_INTERFACE
#undef HAVE_DARWIN_IEEE_INTERFACE
#undef HAVE_DARWIN86_IEEE_INTERFACE])

AH_BOTTOM([/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif])

AH_BOTTOM([/* Substitute gsl functions for missing system functions */

#if !HAVE_DECL_HYPOT
#define hypot gsl_hypot
#endif

#if !HAVE_DECL_LOG1P
#define log1p gsl_log1p
#endif

#if !HAVE_DECL_EXPM1
#define expm1 gsl_expm1
#endif

#if !HAVE_DECL_ACOSH
#define acosh gsl_acosh
#endif

#if !HAVE_DECL_ASINH
#define asinh gsl_asinh
#endif

#if !HAVE_DECL_ATANH
#define atanh gsl_atanh
#endif

#if !HAVE_DECL_LDEXP
#define ldexp gsl_ldexp
#endif

#if !HAVE_DECL_FREXP
#define frexp gsl_frexp
#endif

#if !HAVE_DECL_ISINF
#define isinf gsl_isinf
#endif

#if !HAVE_DECL_FINITE
#define finite gsl_finite
#endif

#if !HAVE_DECL_ISNAN
#define isnan gsl_isnan
#endif])

AH_BOTTOM([#ifdef __GNUC__
#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);
#else
#define DISCARD_POINTER(p) /* ignoring discarded pointer */
#endif])

AH_BOTTOM([#if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#define GSL_RANGE_CHECK 0  /* turn off range checking by default internally */
#endif])

AH_BOTTOM([#define RETURN_IF_NULL(x) if (!x) { return ; }
])

AH_VERBATIM([GSL_DISABLE_DEPRECATED],
[/* Disable deprecated functions and enums while building */
#define GSL_DISABLE_DEPRECATED 1])

dnl
AC_CONFIG_FILES([gsl_version.h gsl.spec gsl/Makefile test/Makefile err/Makefile sys/Makefile utils/Makefile const/Makefile min/Makefile multimin/Makefile ieee-utils/Makefile fft/Makefile specfunc/Makefile dht/Makefile fit/Makefile multifit/Makefile bspline/Makefile statistics/Makefile sum/Makefile roots/Makefile multiroots/Makefile ntuple/Makefile poly/Makefile qrng/Makefile rng/Makefile randist/Makefile siman/Makefile integration/Makefile interpolation/Makefile doc/Makefile block/Makefile vector/Makefile matrix/Makefile histogram/Makefile monte/Makefile ode-initval/Makefile cblas/Makefile blas/Makefile linalg/Makefile eigen/Makefile permutation/Makefile combination/Makefile sort/Makefile complex/Makefile diff/Makefile deriv/Makefile cheb/Makefile cdf/Makefile wavelet/Makefile Makefile])
AC_OUTPUT
