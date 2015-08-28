#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>



#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>
#include <errno.h>
#include "error.h"


#define POUT2(caller, x, a, b)     pout2(caller, #x, x, #a, a, #b, b)
#define POUT3(caller, x, a, b, c)  pout3(caller, #x, x, #a, a, #b, b, #c, c)


static void pout2(const char *caller,
		  const char *sval, int ival, 
		  const char *sval1, int oval1,
		  const char *sval2, int oval2)
{
  if ( ival == oval1 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if ( ival == oval2 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else
    fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}


static void pout3(const char *caller,
		  const char *sval, int ival, 
		  const char *sval1, int oval1,
		  const char *sval2, int oval2,
		  const char *sval3, int oval3)
{
  if ( ival == oval1 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if ( ival == oval2 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else if ( ival == oval3 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval3);
  else
    fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}


void print_pthread_attr(const char *caller, pthread_attr_t *attr)
{
  struct sched_param param;
  int detachstate, policy, inherit, scope, priority;
  size_t stacksize;

  pthread_attr_getdetachstate(attr, &detachstate);
  POUT2(caller, detachstate, PTHREAD_CREATE_JOINABLE, PTHREAD_CREATE_DETACHED);

#if defined(SCHED_FIFO)
  pthread_attr_getschedpolicy(attr, &policy);
  POUT3(caller, policy, SCHED_FIFO, SCHED_RR, SCHED_OTHER);
  pthread_attr_getschedparam(attr, &param);
  priority = param.sched_priority;
  fprintf(stderr, "%-18s :  %-14s = %d\n", caller, "priority", priority);
#endif

#if defined(PTHREAD_INHERIT_SCHED)
  pthread_attr_getinheritsched(attr, &inherit);
  POUT2(caller, inherit, PTHREAD_INHERIT_SCHED, PTHREAD_EXPLICIT_SCHED);
#endif

  pthread_attr_getscope(attr, &scope);
  POUT2(caller, scope, PTHREAD_SCOPE_SYSTEM, PTHREAD_SCOPE_PROCESS);

  pthread_attr_getstacksize(attr, &stacksize);
  fprintf(stderr, "%-18s :  %-14s = %ld\n", caller, "stacksize", (long) stacksize);
}


void print_pthread_mutexattr(const char *caller,  pthread_mutexattr_t *m_attr)
{
  /*
#if defined(_POSIX_THREAD_PRIO_PROTECT) && defined(_POSIX_THREAD_PRIO_INHERIT)
  {
  int protocol;
  pthread_mutexattr_getprotocol(m_attr, &protocol);
  POUT3(caller, protocol, PTHREAD_PRIO_INHERIT, PTHREAD_PRIO_PROTECT, PTHREAD_PRIO_NONE);
  }
#endif
  */
#if defined(PTHREAD_MUTEX_FAST_NP)
  {
  int kind;
  pthread_mutexattr_getkind_np(m_attr, &kind);
  POUT3(caller, kind, PTHREAD_MUTEX_FAST_NP, PTHREAD_MUTEX_RECURSIVE_NP, PTHREAD_MUTEX_ERRORCHECK_NP);
  }
#endif

#if defined(_POSIX_THREAD_PROCESS_SHARED)
  {
  int pshared;
  pthread_mutexattr_getpshared(m_attr, &pshared);
  POUT2(caller, pshared, PTHREAD_PROCESS_SHARED, PTHREAD_PROCESS_PRIVATE);
  }
#endif
}


void print_pthread_condattr(const char *caller, pthread_condattr_t *c_attr)
{
#if defined(_POSIX_THREAD_PROCESS_SHARED)
  {
  int pshared;
  pthread_condattr_getpshared(c_attr, &pshared);
  POUT2(caller, pshared, PTHREAD_PROCESS_SHARED, PTHREAD_PROCESS_PRIVATE);
  }
#endif
}


int PTHREAD_Debug = 0;


void Pthread_debug(int debug)
{
  PTHREAD_Debug = debug;
}


int Pthread_create(const char *caller, pthread_t *th,
		   pthread_attr_t *attr, void * (*start_routine)(void *), void *arg)
{
  int status;

  if ( PTHREAD_Debug ) Message("%s", caller);

  if ( PTHREAD_Debug )
    {
      Message("%s attributes:", caller);
      if ( attr )
	print_pthread_attr(__func__, attr);
      else
	Message("  default attributes");
    }

  status = pthread_create(th, attr, start_routine, arg);

  //if ( PTHREAD_Debug ) Message("-%s (thID = %ld, status = %d)", caller, (long) *th, status);

  return (status);
}


int Pthread_join(const char *caller, pthread_t th, void **thread_return)
{
  int status;

  //  if ( PTHREAD_Debug ) Message("+%s (thID = %ld)", caller, (void *) th);

  status = pthread_join(th, thread_return);

  // if ( PTHREAD_Debug ) Message("-%s (thID = %ld, status = %d)", caller, (void *) th, status);

  return (status);
}


void Pthread_mutex_lock(const char *caller, pthread_mutex_t *mutex)
{
  int status;

  if ( PTHREAD_Debug ) Message("%s (mutex = %p)", caller, (void *) mutex);

  status = pthread_mutex_lock(mutex);
  if ( status != 0 )
    {
      switch (status)
	{
	case EINVAL:
	  Error("The mutex has not been properly initialized!");
	  break;
	case EDEADLK:
	  Error("The mutex is already locked by the calling thread!");
	  break;
	default:
	  Error("Status %d unknown!", status, (void *) mutex);
	  break;
	}
    }
}


void Pthread_mutex_unlock(const char *caller, pthread_mutex_t *mutex)
{
  int status;

  if ( PTHREAD_Debug ) Message("%s (mutex = %p)", caller, (void *) mutex);

  status = pthread_mutex_unlock(mutex);
  if ( status != 0 )
    {
      switch (status)
	{
	case EINVAL:
	  Error("The mutex has not been properly initialized!");
	  break;
	case EPERM:
	  Error("The calling thread does not own the mutex!");
	  break;
	default:
	  Error("Status %d unknown!", status);
	  break;
	}
    }
}


void Pthread_cond_signal(const char *caller, pthread_cond_t *cond)
{
  if ( PTHREAD_Debug ) Message("+%s (cond = %p)", caller, (void *) cond);

  pthread_cond_signal(cond);

  if ( PTHREAD_Debug ) Message("-%s (cond = %p)", caller, (void *) cond);
}


void Pthread_cond_wait(const char *caller, pthread_cond_t *cond, pthread_mutex_t *mutex)
{
  if ( PTHREAD_Debug ) Message("+%s (cond = %p, mutex =  %p)",
			       caller, (void *) cond, (void *) mutex);

  pthread_cond_wait(cond, mutex);

  if ( PTHREAD_Debug ) Message("-%s (cond = %p, mutex = %p)",
			       caller, (void *) cond, (void *) mutex);
}

#endif
