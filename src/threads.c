/** Types for Multi-threading */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2012 - 2014 Genome Research Ltd.                           *
 *                                                                           *
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              *
 *                                                                           *
 *  This file is part of SMALT.                                              *
 *                                                                           *
 *  SMALT is free software: you can redistribute it and/or modify it under   *
 *  the terms of the GNU General Public License as published by the Free     *
 *  Software Foundation, either version 3 of the License, or (at your        *
 *  option) any later version.                                               *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *  General Public License for more details.                                 *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/
#include <stdlib.h>
#include <limits.h>
#include <semaphore.h>
#include <pthread.h>
#include "threads.h"
#ifdef THREADS_DEBUG
#include <string.h>
#include <stdarg.h>
#endif


#ifdef __APPLE__
#define macosx_semaphores
/* use named semaphores for MacOSX Posix threads API 
* unnamed semaphores are not supported */
#endif


enum TRHEAD_CONST {
  BUFFARG_NUM_FAC_DEFAULT  = 4, /**< Number of buffered thread arguments/results
				 * as a multiplicative factor. I.e. for n>0 threads there are
				 * n*BUFFARG_NUM_FAC buffered arguments */
  THREAD_BUFF_NUM = 3,  /**< Number of buffers for arguments that are passed between
			 * threads. */
  THREAD_TASK_NUM = 4,   /**< Number of tasks available */
#ifdef macosx_semaphores
  NAMED_SEMA_PERM = 0644, /**< user: read + write, others/group: read */
#endif
#ifdef THREADS_DEBUG
  DBG_CHARBUFSIZ = 32,    /**< buffer size for debugging output */
#endif
};

enum THREADS_STATUS_FLAGS {
  THRFLG_INIT = 0x01,  /**< Threads were initialised */
  THRFLG_PROC = 0x02,
  THRFLG_INPUT = 0x04,
  THRFLG_OUTPUT = 0x08,
  THRFLG_ARGBUF = 0x10, /**< Argument buffers have been set up */
  THRFLG_SETUP = 0x20,   /**< Threads have been set up */
  THRFLG_STARTED = 0x40
};

enum THREAD_ARGUMENT_FLAGS {
  THRARG_INDEPT = 0x01,
  THRARG_SIGNED = 0x02,
};

enum THREAD_BUFFER_TYPES {
  THRBUFTYP_EMPTY = 0,
  THRBUFTYP_LOADED = 1,
  THRBUFTYP_PROCESSED = 2,
  THRBUFTYP_UNKNOWN = 3,
  THRBUFTYP_NUM = 4
};

typedef uint8_t BOOL_t;
typedef uint8_t THRSTATUSFLG_t;
typedef uint8_t THRARGFLG_t;
typedef uint8_t THRBUFTYP_t; /* one of THREAD_BUFFER_TYPES */
typedef void *(PTHREAD_WRAPPER_FUNC)(void *);

typedef struct _BUFFARG { /**< Arguments passed between threads */
#ifdef THREADS_DEBUG
  short threadno;
  uint64_t readno;
#endif
  short argno;
  void *thisp;
  struct _BUFFARG *nextp;
} BUFFARG;

typedef struct _ARGBUFF { /**< Buffer for arguments to be passed between threads */
  int nThreadsPushing;  /**< Number of threads pushing on this buffer,
			 * 0 signals termination */
  BUFFARG *firstp;
  BUFFARG *lastp;
#ifdef macosx_semaphores
  sem_t *sema;
#else
  sem_t sema; 
  /* semaphore, keeps track of how many arguments are in buffer */
#endif
  pthread_mutex_t mutex; /* lock access to shared data */
  THRBUFTYP_t buftyp;    /* one of THREAD_BUFFER_TYPES */
#ifdef THREADS_DEBUG
  int n_buffarg;
#endif
} ARGBUFF;

typedef struct _THREADTASK {
  THRSTATUSFLG_t status;
  short n_threads;      /**< Number of threads for this task (can be 0) */
  THREAD_INITF *initf;   /**< Initialisation function */
  const void *initargp;  /**< Additional argument to pass to intialisation function */
  THREAD_PROCF *procf;   /**< Processing function */
  THREAD_CLEANF *cleanf; /**< Function for cleaning up */
  THREAD_CHECKF *checkf; /**< Checks buffer before fetching */
  THREAD_CMPF *cmpf;     /**< Compare buffer arguments */
  size_t argsz;         /**< Size of the thread argument */
  uint8_t fromx;        /**< Index of buffer from which arguments are pulled */
  uint8_t tox;          /**< Index of buffer to which argments are pushed */
} THREADTASK;

typedef struct _THREADARG {
  short threadno;   /**< Thread number (< 0 means not run as a separate thread */
  THRARGFLG_t flags; /**< Combination of THREAD_ARGUMENT_FLAGS */
  uint8_t task;     /**< one of THREAD_TASKS */
  void *p;          /**< points to thread-specific data (SmaltMapArgs or SmaltInputArgs) */
  BUFFARG *buflstp;  /**< Start of linked list for internal buffering (e.g. sorted ouptput) */
  ErrMsg *errmsgp;  /**< Thread specific error messages */
  int exit_code; /**< error code with which wrapper exits */
} THREADARG;

static struct _Threads {
  uint8_t status;       /**< One of THREADS_STATUS_FLAGS */
  short n_threads;      /**< Number of threads -1 (0 if run single-threaded) */
  THREADTASK tasks[THREAD_TASK_NUM];
  pthread_t *threadp;   /**< Thread structures, array of size n_threads */
  short n_targ;         /**< Size of array targp */
  THREADARG *targp;     /**< Thread arguments array of size n_targ */
  short n_buffargs;     /**< Number of bufferend arguments */
  BUFFARG *buffargp;    /**< Buffered arguments, array of size n_buffargs */
  void *memp;           /**< Memory allocated for arguments */
  ARGBUFF buff[THREAD_BUFF_NUM];/**< Buffers for arguments that are passed between threads.
  				 * buffp[0]: empty arguments, bufp[1]: loaded arguments.
  				 * buffp[2]: processed arguments. */
} Threads;


#ifdef THREADS_DEBUG
static pthread_mutex_t mutex_dbgout; //= PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef THREADS_DEBUG
static int getBufArgNum(const ARGBUFF *fifop)
{
  int n = 0;
  BUFFARG *hp;
  for (hp = fifop->firstp; NULL != hp; hp = hp->nextp, n++);
  return n;
}

static const char *getTaskTyp(const THREADARG *p)
{
  switch(p->task) {
  case THRTASK_ARGBUF:
    return "argument";
  case THRTASK_INPUT:
    return "input";
  case THRTASK_PROC:
    return "processing";
  case THRTASK_OUTPUT:
    return "output";
  default:
    return "unknown";
  }
}

static const char *getBufTyp(const ARGBUFF *fifop)
{
  switch (fifop->buftyp) {
  case THRBUFTYP_EMPTY:
    return "empty";
  case THRBUFTYP_LOADED:
    return "loaded";
  case THRBUFTYP_PROCESSED:
    return "processed";
  case THRBUFTYP_UNKNOWN:
  default:
    return "unknown";
  }
}
#endif

#ifdef macosx_semaphores
static const char *getBufSemaphorNam(THRBUFTYP_t buftyp)
{
  switch (buftyp) {
  case THRBUFTYP_EMPTY:
    //return "/sem_empty";
    return "SEMAPHORE_BUFF_EMPTY";
  case THRBUFTYP_LOADED:
    //return "/sem_load";
    return "SEMAPHORE_BUFF_LOADED";
  case THRBUFTYP_PROCESSED:
    //return "/sem_proc";
    return "SEMAPHORE_BUFF_PROCESSED";
  case THRBUFTYP_UNKNOWN:
  default:
    //return "/sem_unknown";
    return "SEMAPHORE_BUFF_UNKNOWN";
  }
}
#endif

/*****************************************************************************
 *********************** Methods of private type ARGBUFF *********************
 *****************************************************************************/

static int pushARGBUFF(ARGBUFF *fifop, BUFFARG *argp)
/**< Put argument in buffer. argp == NULL initialises. */
{
  int errcode = ERRCODE_SUCCESS;

  if (argp == NULL) { /* initialisation signal */
    pthread_mutex_lock(&fifop->mutex);
    fifop->firstp = fifop->lastp = NULL;
    fifop->nThreadsPushing = 0;
#ifdef THREADS_DEBUG
    fifop->n_buffarg = 0;
#endif
    pthread_mutex_unlock(&fifop->mutex);
#ifdef THREADS_DEBUG
    threadsPrintDebugMsg(NULL, fifop, "initialised");
#endif
  } else {
    argp->nextp = NULL;
    pthread_mutex_lock(&fifop->mutex);
    if (fifop->firstp == NULL) {
      if (fifop->lastp != NULL)
	errcode = ERRCODE_ASSERT;
      fifop->firstp = fifop->lastp = argp;
#ifdef THREADS_DEBUG
      fifop->n_buffarg = 1;
#endif
    } else {
      if (fifop->lastp == NULL) {
	errcode = ERRCODE_NULLPTR;
      } else {
	fifop->lastp->nextp = argp;
	fifop->lastp = argp;
#ifdef THREADS_DEBUG
	fifop->n_buffarg++;
#endif
      }
    }
    pthread_mutex_unlock(&fifop->mutex);
#ifdef THREADS_DEBUG 
    threadsPrintDebugMsg(argp, fifop,
			 "pushed onto buffer, errcode = %i ...",
			 errcode);
#endif

#ifdef macosx_semaphores 
    sem_post(fifop->sema);
#else
    sem_post(&fifop->sema);
#endif
  }

  return errcode;
}

static int pullARGBUFF(BUFFARG **argp, 
		       ARGBUFF *fifop)
/**< Fetch argument from buffer. Returns termination code 
 *  ERRCODE_PTHRTERMSIG if no argument in buffer and 
 * ifop->nThreadsPushing == 0.
 */
{
  int errcode = ERRCODE_SUCCESS;
#ifdef THREADS_DEBUG
  int na = 0;
#endif

#ifdef THREADS_DEBUG 
  threadsPrintDebugMsg(NULL, fifop, "pulling from buffer, waiting ...");
#endif
#ifdef macosx_semaphores
  sem_wait(fifop->sema);
#else
  sem_wait(&fifop->sema);
#endif
#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, fifop, "pulling from buffer: finished waiting."); 
#endif

  pthread_mutex_lock(&fifop->mutex);
#ifdef THREADS_DEBUG
  na = getBufArgNum(fifop);
#endif
  *argp = fifop->firstp;
  if ((*argp) == NULL) {
    if (fifop->nThreadsPushing <= 0)
      errcode = ERRCODE_PTHRTERMSIG;
    else if (fifop->lastp != NULL)
      errcode = ERRCODE_ASSERT;
  } else {
    if (fifop->firstp == fifop->lastp) {
      if ((*argp)->nextp != NULL)
	errcode = ERRCODE_ASSERT;
      fifop->firstp = fifop->lastp = NULL;
    } else {
      fifop->firstp = (*argp)->nextp;
      (*argp)->nextp = NULL;
    }
#ifdef THREADS_DEBUG
    if (!errcode) fifop->n_buffarg--;
#endif
  }
  pthread_mutex_unlock(&fifop->mutex);
  if (errcode == ERRCODE_PTHRTERMSIG) {
#ifdef macosx_semaphores
    sem_post(fifop->sema);
#else
    sem_post(&fifop->sema); /* release other threads waiting to pull */
#endif
  }
#ifdef THREADS_DEBUG 
  threadsPrintDebugMsg(*argp, fifop,"finished pulling from buffer, "	\
  		       "n_args_prev = %i, errcode = %i ...",
  		       na, errcode);
  if (ERRCODE_PTHRTERMSIG == errcode)
    threadsPrintDebugMsg(argp, fifop, "status: terminated.");
#endif

  return errcode;
}

static int initARGBUFF(ARGBUFF *fifop, THRBUFTYP_t buftyp)
/**< Initialise buffer */
{
  int errcode;
#ifdef macosx_semaphores
  const char * const snam = getBufSemaphorNam(buftyp);
  sem_unlink(snam);
  fifop->sema = sem_open(snam, 
			 O_CREAT, NAMED_SEMA_PERM, 0);
  if (SEM_FAILED == fifop->sema)
    return ERRCODE_SEMOPEN;
#ifdef THREADS_DEBUG
  fifop->buftyp = buftyp;
  threadsPrintDebugMsg(NULL, fifop, "opening named semaphore %s, %p", 
		       snam, fifop->sema);
#endif

#else
  sem_init(&fifop->sema, 0, 0);
#endif
  fifop->buftyp = buftyp;
  pthread_mutex_init(&fifop->mutex, NULL);
  errcode = pushARGBUFF(fifop, NULL);
  return errcode;
}

static void cleanupARGBUFF(ARGBUFF *fifop)
{
  pthread_mutex_destroy(&fifop->mutex);
#ifdef macosx_semaphores
  sem_close(fifop->sema);
#endif
}

static int signOnARGBUFF(ARGBUFF *fifop)
/**< Increment the counter for threads pushing onto the buffer and
 * return the number */
{
  int n;
  pthread_mutex_lock(&fifop->mutex);
  n = ++(fifop->nThreadsPushing);
  pthread_mutex_unlock(&fifop->mutex);
#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, fifop, "signed on ...");
#endif
  return n;
}

static int signOffARGBUFF(ARGBUFF *fifop)
/**< Decrement the counter for threads pushing onto the buffer.
 * Return the number of threads still pushing onto the buffer.
 * If there is no thread pushing, wake threads waiting to pull from the buffer.
 */
{
  int n;
  pthread_mutex_lock(&fifop->mutex);
  n = --(fifop->nThreadsPushing);
  pthread_mutex_unlock(&fifop->mutex);
  if (n <= 0) {
#ifdef macosx_semaphores
    sem_post(fifop->sema);
#else
    sem_post(&fifop->sema);
#endif
  }
  
#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, fifop, "signed off ...", getBufTyp(fifop), n);
#endif

  return n;
}

#ifdef THREADS_SUPERFLUOUS
static void termARGBUFF(ARGBUFF *fifop)
/**< Set termination signal for buffer */
{
#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, fifop, "termARGBUFF() called ..."); 
#endif
  pthread_mutex_lock(&fifop->mutex);
  fifop->nThreadsPushing = 0;
  pthread_mutex_unlock(&fifop->mutex);
#ifdef macosx_semaphores
    sem_post(fifop->sema);
#else
    sem_post(&fifop->sema);
#endif
}
#endif

/****************************************************************************
 ******************** Methods of Private Type THREADARG *********************
 ****************************************************************************/

static void pushTHREADARGInternalBuffer(THREADARG *thargp, BUFFARG *argp, 
				       THREAD_CMPF *cmpf)
{
  BUFFARG *hp = thargp->buflstp;
#ifdef THREADS_DEBUG
  BUFFARG *p;
  int i = 0;
  if (NULL != cmpf) {
    pthread_mutex_lock(&mutex_dbgout);
  
    fprintf(stderr, "THREADS_DEBUG: +++++ pushTHREADARGInternalBuffer: "\
	    "threadno = %hi, readno = %llu +++++\n",
	   thargp->threadno, (unsigned long long) argp->readno); 
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      fprintf(stderr, 
	      "THREADS_DEBUG: pushTHREADARGInternalBuffer: [%i] %llu\n",
	      i, (unsigned long long) p->readno);
    pthread_mutex_unlock(&mutex_dbgout);
  }
#endif
  if (NULL == cmpf || NULL == hp || 
      cmpf(argp->thisp, hp->thisp) <= 0) {
    argp->nextp = hp;
    thargp->buflstp = argp;
  } else {
    BUFFARG *lp = hp->nextp;
    for(; (lp) && cmpf(argp->thisp, lp->thisp) > 0; lp = lp->nextp)
      hp = lp;
    argp->nextp = lp;
    hp->nextp = argp;
  }
#ifdef THREADS_DEBUG
  if (NULL != cmpf) {
    pthread_mutex_lock(&mutex_dbgout);
    fprintf(stderr,
	    "THREADS_DEBUG: + pushTHREADARGInternalBuffer: " \
	    "on exit: threadno = %hi +\n", thargp->threadno);
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      fprintf(stderr,
	      "THREADS_DEBUG: pushTHREADARGInternalBuffer: readno[%i] = %llu\n",
	     i, (unsigned long long) p->readno);
    fprintf(stderr,"THREADS_DEBUG: ++ pushTHREADARGInternalBuffer: exit ++\n");
    pthread_mutex_unlock(&mutex_dbgout);
  }
#endif
}

static BUFFARG *pullTHREADARGInternalBuffer(THREADARG *thargp, THREAD_CHECKF *checkf, void *tdatap)
{
  BUFFARG *argp = thargp->buflstp;
#ifdef THREADS_DEBUG
  int chkrv = 0;
  BUFFARG *p;
  int i = 0;
  if (NULL != checkf) {
    pthread_mutex_lock(&mutex_dbgout);
    fprintf(stderr,
	   "THREADS_DEBUG: ----- pullTHREADARGInternalBuffer "\
	   "threadno = %hi-----\n", thargp->threadno); 
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      fprintf(stderr,
	      "THREADS_DEBUG: pullTHREADARGInternalBuffer: [%i] %llu\n",
	     i, (unsigned long long) p->readno);
    pthread_mutex_unlock(&mutex_dbgout);
  }
#endif

  if (argp != NULL) {

    if (NULL == checkf ||
#ifdef THREADS_DEBUG
	(chkrv = checkf(tdatap, argp->thisp))) {
#else
      (checkf(tdatap, argp->thisp))) {
#endif
      thargp->buflstp = argp->nextp;
      argp->nextp = NULL;
    } else {
      argp = NULL;
    }
  }

#ifdef THREADS_DEBUG
  if (NULL != checkf) {
    pthread_mutex_lock(&mutex_dbgout);
    fprintf(stderr, "THREADS_DEBUG: pullTHREADARGInternalBuffer: checkf = %i\n",
	   chkrv);
    fprintf(stderr, 
	    "THREADS_DEBUG: - pullTHREADARGInternalBuffer: "	\
	    "on exit: threadno = %hi -\n", thargp->threadno);
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      fprintf(stderr,
	      "THREADS_DEBUG: pullTHREADARGInternalBuffer: readno[%i] = %llu\n",
	     i, (unsigned long long) p->readno);
    fprintf(stderr,
	    "THREADS_DEBUG: -- pullTHREADARGInternalBuffer: exit --\n");
    pthread_mutex_unlock(&mutex_dbgout);
  }
#endif

  return argp;
}
/****************************************************************************
 *************** PTHREAD_WRAPPER_FUNC function run by threads ***************
 ****************************************************************************/

static void *tprocf(void *p) 
{
  int errcode = ERRCODE_SUCCESS;
  BUFFARG *argp;
  THREADARG *thargp = (THREADARG *) p;
  THREADTASK *taskp = Threads.tasks + thargp->task;
  ARGBUFF *argbf_fromp = Threads.buff + taskp->fromx;
  ARGBUFF *argbf_top = Threads.buff + taskp->tox;

#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, NULL, "thead[%hi]: task '%s' started",
  		       thargp->threadno, getTaskTyp(thargp));
#endif

  if ( !(thargp->flags & THRARG_SIGNED)) {
#ifdef THREADS_DEBUG
    threadsPrintDebugMsg(NULL, Threads.buff + taskp->tox, 
			 "thead[%hi]:tprocf(): task '%s': sign on", 
			 thargp->threadno, getTaskTyp(thargp));
#endif
    signOnARGBUFF(argbf_top);
    thargp->flags |= THRARG_SIGNED;
  }

  do {

    if ((errcode = pullARGBUFF(&argp, argbf_fromp)))
      break;

#ifdef THREADS_DEBUG
    threadsPrintDebugMsg(argp, argbf_fromp, 
			 "thead[%hi]:tprocf(): pulled argument\n",
    			 thargp->threadno);
#endif
    if (argp == NULL) {
      errcode = ERRCODE_ASSERT;
      break;
    } 
    
#ifdef THREADS_DEBUG
    threadsPrintDebugMsg(argp, NULL, 
			 "thead[%hi]:tprocf('%s', %hu): "\
			 "push to internal buffer cmpf %c= NULL",
    			 thargp->threadno, getTaskTyp(thargp),
    			 (unsigned short) thargp->task,
    			 (NULL == taskp->cmpf)? '=':'!');
#endif

    pushTHREADARGInternalBuffer(thargp, argp, taskp->cmpf);

    while ((ERRCODE_SUCCESS == errcode) && 
	   (argp = pullTHREADARGInternalBuffer(thargp, taskp->checkf, thargp->p))) {
      errcode = taskp->procf(thargp->errmsgp, 
#ifdef THREADS_DEBUG
			     &argp->readno,
#endif
			     thargp->p, argp->thisp);

    /* push to loaded arguments */
#ifdef THREADS_DEBUG
      threadsPrintDebugMsg(argp, argbf_top,
			   "thread[%hi]:tprocf push read onto buffer %hi",
      			   thargp->threadno, taskp->tox);
#endif
      if (!(errcode))
	errcode = pushARGBUFF(argbf_top, argp);
    }
 
  } while ((ERRCODE_SUCCESS == errcode)
	   && thargp->threadno >= 0);

  if ((errcode)) { 
#ifdef THREADS_DEBUG
    threadsPrintDebugMsg(NULL, argbf_top, 
			 "thead[%hi]:sign off from buffer %hi (errcode = %i).",
    			 thargp->threadno, taskp->tox, errcode);
#endif

    signOffARGBUFF(argbf_top);
    thargp->flags &= ~THRARG_SIGNED;
  }
#ifdef THREADS_DEBUG
  threadsPrintDebugMsg(NULL, argbf_top,
		      "thead[%hi]:tprocf(): task '%s' finished.\n",
  		       thargp->threadno, getTaskTyp(thargp));
#endif
  thargp->exit_code = errcode;

  return p;
}

/****************************************************************************
 ***************************** Public Methods *******************************
 ****************************************************************************/
#ifdef THREADS_DEBUG

int threadsPrintDebugMsg(void *argp, void *fifop, const char *format, ...)
{
  int nchar = 0;
  int nt = 0;
  int na = 0;
  int ns = 0;
  short threadno = 0;
  int argno = 0;
  int readno = 0;
  char chrbf[DBG_CHARBUFSIZ];

  if (argp) {
    BUFFARG * const p = (BUFFARG *) argp;
    threadno = p->threadno;
    argno = p->argno;
    readno = p->readno;
  }
  if (fifop) {
    ARGBUFF * const p =  (ARGBUFF *) fifop;
    pthread_mutex_lock(&p->mutex);
    nt = p->nThreadsPushing;
    na = getBufArgNum(p);
#ifdef macosx_semaphores
    ns = p->n_buffarg;
#else
    sem_getvalue(&p->sema, &ns);
#endif
    strncpy(chrbf, getBufTyp(p), DBG_CHARBUFSIZ);
    chrbf[DBG_CHARBUFSIZ-1]='\0';
    pthread_mutex_unlock(&p->mutex);
  } else {
    chrbf[0] = '\0';
  }

  pthread_mutex_lock(&mutex_dbgout);
  fprintf(stderr, "#THREADS_DEBUG: ");
  if (argp) 
    fprintf(stderr, "thread[%hi] (arg %i, read %i) ", 
	    threadno, argno, readno);
  if (fifop) 
    fprintf(stderr, "'%s' buffer "\
	    "(nThreadsPushing = %i, n_args = %i, sema = %i):\n", 
	    chrbf, nt, na, ns);
  
  if (format != NULL) {
    va_list ap;
    va_start(ap, format);
    nchar = vfprintf(stderr, format, ap);
    va_end(ap);
  }
  fprintf(stderr, "\n");

  pthread_mutex_unlock(&mutex_dbgout);

  return nchar;
}
#endif

int threadsInit(void)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  Threads.status = 0;
  Threads.threadp = NULL;
  Threads.targp = NULL;
  Threads.buffargp = NULL;
  for (i=0; i<THREAD_BUFF_NUM && ERRCODE_SUCCESS == errcode; i++) {
    THRBUFTYP_t buftyp = THRBUFTYP_UNKNOWN;
    if (i == 0)
      buftyp = THRBUFTYP_EMPTY;
    else if (i == 1)
      buftyp = THRBUFTYP_LOADED;
    else if (i == 2)
      buftyp = THRBUFTYP_PROCESSED;
    else
      buftyp = THRBUFTYP_UNKNOWN;     
    errcode = initARGBUFF(Threads.buff + i, buftyp);
  }
  Threads.status |= THRFLG_INIT;
  return errcode;
}

int threadsSetTask(uint8_t task_typ, 
		   short n_threads, 
		   THREAD_INITF *initf, 
		   const void *initargp,
		   THREAD_PROCF *procf, 
		   THREAD_CLEANF *cleanf,
		   THREAD_CHECKF *checkf,
		   THREAD_CMPF *cmpf,
		   size_t argsz)
{
  int errcode = ERRCODE_SUCCESS;
  THREADTASK *taskp = NULL;

  if (task_typ >= THREAD_TASK_NUM)
    return ERRCODE_FAILURE;

  if (!(Threads.status & THRFLG_INIT))
    return ERRCODE_ASSERT;

  taskp = Threads.tasks + task_typ;
  taskp->n_threads = n_threads;
  taskp->initf = initf;
  taskp->initargp = initargp;
  taskp->procf = procf;
  taskp->cleanf = cleanf;
  taskp->checkf = checkf;
  taskp->cmpf = cmpf;
  taskp->argsz = argsz;
  taskp->status = THRFLG_INIT;
  taskp->fromx = taskp->tox = 0;

  switch (task_typ) {
  case THRTASK_INPUT:
    taskp->fromx = THRBUFTYP_EMPTY;
    taskp->tox = THRBUFTYP_LOADED;
    Threads.status |= THRFLG_INPUT;
    break;
  case THRTASK_ARGBUF:
    taskp->n_threads = 0;
    taskp->procf = NULL;
    Threads.status |= THRFLG_ARGBUF;
    break;
  case THRTASK_PROC:
    taskp->fromx = THRBUFTYP_LOADED;
    taskp->tox = THRBUFTYP_PROCESSED;
    Threads.status |= THRFLG_PROC;
    break;
  case THRTASK_OUTPUT:
    taskp->fromx = THRBUFTYP_PROCESSED;
    taskp->tox = THRBUFTYP_EMPTY;
    Threads.status |= THRFLG_OUTPUT;
#ifdef THREADS_DEBUG
    fprintf(stderr,
	    "#THREADS_DEBUG: threadsSetTask [%hu]: cmpf %c= NULL, checkf %c= NULL\n",
	   (unsigned short) task_typ, 
	   (NULL ==  taskp->cmpf)? '=':'!', 
	   (NULL ==  taskp->checkf)? '=':'!');
#endif
    break;
  default:
      errcode = ERRCODE_FAILURE;
    break;
  }

  return errcode;
}

int threadsSetUp(int n_buffarg_factor)
{
  int errcode = ERRCODE_SUCCESS;
  uint8_t testflg = THRFLG_INIT | THRFLG_PROC | THRFLG_INPUT | THRFLG_OUTPUT | THRFLG_ARGBUF;
  uint8_t nta;
  short i, nth, ntr, n_threads, n_targ;
  size_t memsz;
  char *hp;

  if (n_buffarg_factor <= 0)
    n_buffarg_factor = BUFFARG_NUM_FAC_DEFAULT;
  if ((Threads.status & testflg) != testflg || 
      ((testflg & THRFLG_SETUP)))
    return ERRCODE_ASSERT;

  Threads.n_threads = 0;
  Threads.n_targ = 0;
  memsz = 0;
  for (nta=0; nta<THREAD_TASK_NUM; nta++) {
    if (nta != THRTASK_ARGBUF) {
      THREADTASK *taskp = Threads.tasks + nta;
      if (taskp->n_threads < 1) {
	Threads.n_targ += 1;
	memsz += taskp->argsz;	
      } else {
	if (Threads.n_threads + taskp->n_threads > SHRT_MAX)
	  return ERRCODE_OVERFLOW;
	Threads.n_threads = (short) (Threads.n_threads + taskp->n_threads);
	Threads.n_targ = (short) (Threads.n_targ + taskp->n_threads);
	memsz += taskp->argsz * taskp->n_threads;
      }
    }
  }

  Threads.n_buffargs = (short) ((Threads.n_threads > 0)? Threads.n_threads * n_buffarg_factor: 1);
  memsz += Threads.n_buffargs * Threads.tasks[THRTASK_ARGBUF].argsz;

#ifdef THREADS_DEBUG
  fprintf(stderr, 
	  "THREADS_DEBUG::threadsSetup(): initialising %hi threads with %hi arguments and %hi buffers\n",
	  Threads.n_threads, Threads.n_targ, Threads.n_buffargs);
#endif

  if (Threads.n_threads > 0) {
    ECALLOCP(Threads.n_threads, Threads.threadp);
    if (NULL == Threads.threadp)
      errcode = ERRCODE_NOMEM;
  } else {
    Threads.threadp = NULL;
  }
    
  ECALLOCP(Threads.n_targ, Threads.targp);  
  ECALLOCP(Threads.n_buffargs, Threads.buffargp);
  Threads.memp = EMALLOC0(memsz);
    
  if (NULL == Threads.targp || 
      NULL == Threads.buffargp ||
      NULL == Threads.memp)
    errcode = ERRCODE_NOMEM;

  n_threads = 0;
  n_targ = 0;
  hp = Threads.memp;
  for (nta=0, nth=0, ntr=0; nta<THREAD_TASK_NUM && !(errcode); nta++) {
    if (nta != THRTASK_ARGBUF) {
      THREADTASK *taskp = Threads.tasks + nta;
      if (taskp->n_threads < 1) {
	n_targ++;
      } else {
	if (n_threads + taskp->n_threads > SHRT_MAX)
	  return ERRCODE_OVERFLOW;
	n_threads = (short) (n_threads + taskp->n_threads);
	n_targ = (short) (n_targ + taskp->n_threads);
      }
      for (; ntr<n_targ && !(errcode); ntr++) {
	THREADARG *targp = Threads.targp + ntr;
	targp->task = nta;
	targp->p = hp;
	targp->buflstp = NULL;
	hp += taskp->argsz;
	if (taskp->n_threads < 1) {
	  targp->threadno = -1;
	  targp->flags = 0;
	} else {
	  targp->threadno = nth++;
	  targp->flags = THRARG_INDEPT;
	}
	if (NULL == ERRMSG_CREATE(targp->errmsgp))
	  errcode = ERRCODE_NOMEM;
	else
	  errcode = (*Threads.tasks[nta].initf)(targp->p, taskp->initargp, targp->threadno);
      }
    }
  }
  if (!(errcode) && 
      (Threads.n_threads != n_threads ||
       Threads.n_targ != n_targ))
    errcode = ERRCODE_ASSERT;

  for (i=0;i<Threads.n_buffargs && (!errcode); i++) {
    BUFFARG *p = Threads.buffargp + i;
    p->thisp = hp;
    hp += Threads.tasks[THRTASK_ARGBUF].argsz;
    p->argno = i;
    p->nextp = NULL;
#ifdef THREADS_DEBUG
    p->threadno = 0;
    p->readno = 0;
#endif
    errcode = (*Threads.tasks[THRTASK_ARGBUF].initf)(p->thisp, 
						      Threads.tasks[THRTASK_ARGBUF].initargp, 
						      i);
    if (ERRCODE_SUCCESS == errcode) {
      errcode = pushARGBUFF(&Threads.buff[THRBUFTYP_EMPTY], p);
    }
  }

  if (ERRCODE_SUCCESS == errcode)
    Threads.status |= THRFLG_SETUP;

  return errcode;
}

void threadsCleanup(void)
{
  if (Threads.status & THRFLG_SETUP) {
    short i;
    for (i=0;i<Threads.n_buffargs; i++) {
      BUFFARG *bargp = Threads.buffargp + i;
      (*Threads.tasks[THRTASK_ARGBUF].cleanf)(NULL, bargp->thisp);
    }
    for (i=0; i<Threads.n_targ; i++) {
      THREADARG *targp = Threads.targp + i;
      (*Threads.tasks[targp->task].cleanf)(targp->errmsgp, targp->p);
      ERRMSG_END(targp->errmsgp)
    }

    for (i=0;i<THREAD_BUFF_NUM; i++) {
      cleanupARGBUFF(Threads.buff + i);
    }
#ifdef maxosx_semaphores
    for (i = THRBUFTYP_NUM-1; i>=0; i--)
      sem_unlink(getBufSemaphorNam(i));
#endif
    free(Threads.memp);
    free(Threads.buffargp);
    free(Threads.targp);
    free(Threads.threadp);
  }
  Threads.status = 0;
}

int threadsStart(void)
{
  int rv = ERRCODE_SUCCESS;
  short ntr;

  if (!(Threads.status & THRFLG_SETUP))
    return ERRCODE_ASSERT;

  /* start up thread for input */ 
  for (ntr=0; ntr < Threads.n_targ; ntr++) {
    THREADARG *targp = Threads.targp + ntr;
    if (targp->threadno < 0)
      continue;
    if ((rv = pthread_create(Threads.threadp + targp->threadno, NULL, 
			     tprocf, targp)))
      ERRMSGNO(targp->errmsgp, ERRCODE_PTHREAD);
  }
  if (!rv) 
    Threads.status |= THRFLG_STARTED;

  return rv;
}

void threadsStop(void)
{
  int nth;
    /* wait for threads to finish */
  for (nth=0; nth < Threads.n_threads; nth++) {
#ifdef THREADS_DEBUG
    pthread_mutex_lock(&mutex_dbgout);
    fprintf(stderr,
	    "# THREAD_DEBUG: main_thread: waiting for thread %hi to finish ...\n", nth);
    pthread_mutex_unlock(&mutex_dbgout);
#endif
    pthread_join(Threads.threadp[nth], NULL);
  }

  return;
}

int threadsRun(void)
{
  int errcode = ERRCODE_SUCCESS;
  BOOL_t has_nonthread = 1;
#ifdef THREADS_DEBUG
  pthread_mutex_init(&mutex_dbgout, NULL);
#endif
  threadsStart();

  
  while (ERRCODE_SUCCESS == errcode && (has_nonthread)) {
    short ntr;
    has_nonthread = 0;
    for(ntr=0; ntr<Threads.n_targ && (ERRCODE_SUCCESS == errcode); ntr++) {
      THREADARG *targp = Threads.targp + ntr;
      if (targp->threadno < 0) {
	has_nonthread = 1;
	tprocf(targp);
	errcode = targp->exit_code;
      }
    }
  }
   
  threadsStop();
#ifdef THREADS_DEBUG
  pthread_mutex_destroy(&mutex_dbgout);
#endif

  return (ERRCODE_PTHRTERMSIG == errcode)? ERRCODE_SUCCESS: errcode;
}

void *threadsGetMem(uint8_t task_typ)
{
  void *p = NULL;

  if (Threads.status & THRFLG_SETUP) {
    if (task_typ == THRTASK_ARGBUF) {
      p = Threads.buffargp[0].thisp;
    } else {
      short ntr;
      for (ntr=0; ntr < Threads.n_targ && p == NULL; ntr++) {
	if (Threads.targp[ntr].task == task_typ) 
	  p = Threads.targp[ntr].p;
      }
    }
  }
  
  return p;
}
