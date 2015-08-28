/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#if defined(HAVE_LIBPTHREAD)
#include <limits.h>
#include <pthread.h>
#include "pthread_debug.h"
#endif


void print_pthread_info()
{
#if defined(HAVE_LIBPTHREAD)
  pthread_attr_t attr;
  pthread_mutexattr_t m_attr;
  pthread_condattr_t c_attr;

#if defined(PTHREAD_KEYS_MAX)
  fprintf(stderr, "PTHREAD_KEYS_MAX    = %d\n", PTHREAD_KEYS_MAX);
#endif
#if defined(PTHREAD_DESTRUCTOR_ITERATIONS)
  fprintf(stderr, "PTHREAD_DESTRUCTOR_ITERATIONS = %d\n", PTHREAD_DESTRUCTOR_ITERATIONS);
#endif
#if defined(PTHREAD_THREADS_MAX)
  fprintf(stderr, "PTHREAD_THREADS_MAX = %d\n", PTHREAD_THREADS_MAX);
#endif
#if defined(PTHREAD_STACK_MIN)
  fprintf(stderr, "PTHREAD_STACK_MIN   = %d\n", PTHREAD_STACK_MIN);
#endif

  fprintf(stderr, "\n");

  pthread_attr_init(&attr);
  print_pthread_attr("Default pthread attr", &attr);
  pthread_attr_destroy(&attr);

  pthread_mutexattr_init(&m_attr);
  print_pthread_mutexattr("Default pthread mutexattr", &m_attr);
  pthread_mutexattr_destroy(&m_attr);

  pthread_condattr_init(&c_attr);
  print_pthread_condattr("Default pthread condattr ", &c_attr);
  pthread_condattr_destroy(&c_attr);

  fprintf(stderr, "\n");
#endif
}
