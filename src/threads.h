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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef THREADS_H
#define THREADS_H

  //#define THREADS_DEBUG

#include <stdint.h>
#include "elib.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/
  
  enum THREAD_TASKS {
    THRTASK_ARGBUF = 0,
    THRTASK_INPUT = 1,
    THRTASK_PROC = 2,
    THRTASK_OUTPUT = 3,
  };

  /****************************************************************************
   ************************** Transparent Types *******************************
   ****************************************************************************/

  typedef void THREADSPEC_t;
  typedef void THREADUNIV_t;
  typedef void THREADEXCH_t;

  /****************************************************************************
   ***************************** Opaque Types *********************************
   ****************************************************************************/

  /****************************************************************************
   *************************** Function Types *********************************
   ****************************************************************************/

  typedef int (THREAD_PROCF) (ErrMsg *,
#ifdef THREADS_DEBUG
			      uint64_t *,
#endif 
			      void *, void *);
  typedef int (THREAD_INITF) (void *, const void *, short);
  typedef int (THREAD_CLEANF)(ErrMsg *, void *);
  typedef int (THREAD_CHECKF)(const void *, const void *);
  typedef int (THREAD_CMPF)(const void *, const void *);

  /****************************************************************************
   ************************ Methods Aiding Debugging **************************
   ****************************************************************************/

#ifdef THREADS_DEBUG
  int threadsPrintDebugMsg(void *, void *, const char *format, ...);
  /**< Prints a message to stderr
   * first 2 arguments used internally.
   */
#endif

  /****************************************************************************
   ********************** Methods of Type ThreadStore *************************
   ****************************************************************************/

  int threadsInit(void);
  /**< Initializer. Has to be called before anything else.
   */
  
  int threadsSetTask(uint8_t task_typ, short n_threads, 
		     THREAD_INITF *initf,
		     const void *initargp,
		     THREAD_PROCF *procf, 
		     THREAD_CLEANF *cleanf,
		     THREAD_CHECKF *checkf,
		     THREAD_CMPF *cmpf,
		     size_t argsz);
  /**< Set functions and parameters for one of the THREAD_TASK tasks
   * \param task_typ One of THREAD_TASK.
   * \param n_threads Number of threads devoted to task,
   * \param initf Initialisation function for task,
   * \param initargp Additional argument to pass to initialisation function. 
   * \param initf Processing function for task,
   * \param cleanf Function that frees allocated memory 
   */

  

  int threadsSetUp(int n_buff_fact);
  /**< Allocate memory and initialise buffers etc.
   * Call this function after threadsSetTask has been called for each of the
   * THREAD_TASKS tasks.
   * \param n_buff_fact Factor determining the buffer size for arguments passed between
   *        threads. The number of elements buffered n_elem = n_threads * n_buff_fact.
   *        n_buff_fact can be <=0 in which case defaults are used.
   */

  void threadsCleanup(void);
  /**< Free all allocated memory, reset buffers etc.
   */

  int threadsStart(void);
  /**< Start all threads.
   */

  void threadsStop(void);
  /**< Stop all threads. Don't clean up. 
   */

  int threadsRun(void);
  /**< Run tasks that were set up serially or in parallel using 
   * multiple threads.
   */

  void *threadsGetMem(uint8_t task_typ);
  /**< Return a pointer to memory associated with a task.
   * \param task_typ One of THREAD_TASKS.
   */

#endif
#ifdef __cplusplus
}
#endif
