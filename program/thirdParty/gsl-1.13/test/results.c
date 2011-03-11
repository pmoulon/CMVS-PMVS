/* err/test_results.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_machine.h>

#if HAVE_VPRINTF
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#endif

#include <gsl/gsl_test.h>

static unsigned int tests = 0;
static unsigned int passed = 0;
static unsigned int failed = 0;

static unsigned int verbose = 0;

static void
initialise (void)
{
  const char * p = getenv("GSL_TEST_VERBOSE");

  /* 0 = show failures only (we always want to see these) */
  /* 1 = show passes and failures */

  if (p == 0)  /* environment variable is not set */
    return ;

  if (*p == '\0') /* environment variable is empty */
    return ;

  verbose = strtoul (p, 0, 0);  

  return;
}

static void 
update (int s)
{
  tests++;

  if (s == 0) 
    {
      passed++;
    }
  else
    {
      failed++;
    }
}

void
gsl_test (int status, const char *test_description,...)
{
  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status && !verbose)
        printf(" [%u]", tests);

      printf("\n");
      fflush (stdout);
    }
}


void
gsl_test_rel (double result, double expected, double relative_error,
              const char *test_description,...)
{
  int status ;

  if (!tests) initialise();

  /* Check for NaN vs inf vs number */

  if (gsl_isnan(result) || gsl_isnan(expected)) 
    {
      status = gsl_isnan(result) != gsl_isnan(expected); 
    }
  else if (gsl_isinf(result) || gsl_isinf(expected)) 
    {
      status = gsl_isinf(result) != gsl_isinf(expected); 
    }
  else if ((expected > 0 && expected < GSL_DBL_MIN)
           || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else if (expected != 0 ) 
    {
      status = (fabs(result-expected)/fabs(expected) > relative_error) ;
    }
  else
    {
      status = (fabs(result) > relative_error) ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
gsl_test_abs (double result, double expected, double absolute_error,
              const char *test_description,...)
{
  int status ;

  if (!tests) initialise();

  /* Check for NaN vs inf vs number */

  if (gsl_isnan(result) || gsl_isnan(expected)) 
    {
      status = gsl_isnan(result) != gsl_isnan(expected); 
    }
  else if (gsl_isinf(result) || gsl_isinf(expected)) 
    {
      status = gsl_isinf(result) != gsl_isinf(expected); 
    }
  else if ((expected > 0 && expected < GSL_DBL_MIN)
           || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else 
    {
      status = fabs(result-expected) > absolute_error ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}


void
gsl_test_factor (double result, double expected, double factor,
                 const char *test_description,...)
{
  int status;

  if (!tests) initialise();
  
  if ((expected > 0 && expected < GSL_DBL_MIN)
      || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else if (result == expected) 
    {
      status = 0;
    }
  else if (expected == 0.0) 
    {
      status = (result > expected || result < expected);
    }
  else
    {
      double u = result / expected; 
      status = (u > factor || u < 1.0 / factor) ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
gsl_test_int (int result, int expected, const char *test_description,...)
{
  int status = (result != expected) ;

  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }
      else 
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n");
      fflush (stdout);
    }
}

void
gsl_test_str (const char * result, const char * expected, 
              const char *test_description,...)
{
  int status = strcmp(result,expected) ;

  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status)
        {
          printf(" (%s observed vs %s expected)", result, expected) ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n");
      fflush (stdout);
    }
}

void
gsl_test_verbose (int v)
{
  verbose = v;
}

int
gsl_test_summary (void)
{
  if (verbose && 0)             /* FIXME: turned it off, this annoys me */
    printf ("%d tests, passed %d, failed %d.\n", tests, passed, failed);

  if (failed != 0)
    {
      return EXIT_FAILURE;
    }

  if (tests != passed + failed)
    {
      if (verbose)
        printf ("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
                tests, passed, failed);
      return EXIT_FAILURE;
    }

  if (passed == tests)
    {
      if (!verbose)         /* display a summary of passed tests */
        printf ("Completed [%d/%d]\n", passed, tests);

      return EXIT_SUCCESS;
    }

  return EXIT_FAILURE;
}
