/**
 * Copyright (C) 2015 Roman Hiestand
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial
 * portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
 * LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "rwmutex.h"
#include <limits.h>

RWMutex::RWMutex()
	: initialized(false), 
	nSharedAccessCount(0),
	nExclusiveAccessCount(0),
	nCompletedSharedAccessCount(0)
{
}

RWMutex::RWMutex(const RWMutex &o)
{
	copy(o);
}

RWMutex::~RWMutex()
{
	destroy();
}

RWMutex &RWMutex::copy(const RWMutex &o)
{
	initialized = o.initialized;
	nSharedAccessCount = o.nSharedAccessCount;
	nExclusiveAccessCount = o.nExclusiveAccessCount;
	nCompletedSharedAccessCount = o.nCompletedSharedAccessCount;
	if(initialized)
	{
		mtxExclusiveAccess = o.mtxExclusiveAccess;
		mtxSharedAccessCompleted = o.mtxSharedAccessCompleted;
		cndSharedAccessCompleted = o.cndSharedAccessCompleted;
	}
	return *this;
}

RWMutex & RWMutex::operator=(const RWMutex & o)
{
	return copy(o);
}

bool RWMutex::init()
{
	if(!initialized)
	{
		mtx_init(&mtxExclusiveAccess, mtx_plain | mtx_recursive);
		mtx_init(&mtxSharedAccessCompleted, mtx_plain | mtx_recursive);
		cnd_init(&cndSharedAccessCompleted);
		initialized = true;
	}
	return true;
}

bool RWMutex::destroy()
{
	if(initialized)
	{
		if(mtx_lock(&mtxExclusiveAccess) != thrd_success)
			return false;

		if(mtx_lock(&mtxSharedAccessCompleted) != thrd_success)
		{
			mtx_unlock(&mtxExclusiveAccess);
			return false;
		}
		if(nExclusiveAccessCount > 0
			|| nSharedAccessCount > nCompletedSharedAccessCount)
		{
			// Busy
			int j = 27;
		}

		mtx_unlock(&mtxSharedAccessCompleted);
		mtx_unlock(&mtxExclusiveAccess);

		cnd_destroy(&cndSharedAccessCompleted);
		mtx_destroy(&mtxSharedAccessCompleted);
		mtx_destroy(&mtxExclusiveAccess);

		initialized = false;
	}
	return true;
}

bool RWMutex::rdlock()
{
	if(mtx_lock(&mtxExclusiveAccess) != thrd_success)
		return false;

//	nSharedAccessCount++;
	if(++nSharedAccessCount == INT_MAX)
	{
		if(mtx_lock(&mtxSharedAccessCompleted) != thrd_success)
		{
			mtx_unlock(&mtxExclusiveAccess);
			return false;
		}
		nSharedAccessCount -= nCompletedSharedAccessCount;
		nCompletedSharedAccessCount = 0;
		if(mtx_unlock(&mtxSharedAccessCompleted) != thrd_success)
		{
			mtx_unlock(&mtxExclusiveAccess);
			return false;
		}
	}

	return (mtx_unlock(&mtxExclusiveAccess) == thrd_success);
}

bool RWMutex::wrlock()
{
	bool result = true;
	if(mtx_lock(&mtxExclusiveAccess) != thrd_success)
		return false;

	if(mtx_lock(&mtxSharedAccessCompleted) != thrd_success)
	{
		mtx_unlock(&mtxExclusiveAccess);
		return false;
	}

	if(nExclusiveAccessCount == 0)
	{
		if(nCompletedSharedAccessCount > 0)
		{
			nSharedAccessCount -= nCompletedSharedAccessCount;
			nCompletedSharedAccessCount = 0;
		}
		if(nSharedAccessCount > 0)
		{
			nCompletedSharedAccessCount = -nSharedAccessCount;

			// pthread_cleanup_push (ptw32_rwlock_cancelwrwait, (void *) rwl);

			do
			{
				result = (cnd_wait(&cndSharedAccessCompleted, &mtxSharedAccessCompleted) == thrd_success);
			}
			while(result && (nCompletedSharedAccessCount < 0));

			// pthread_cleanup_pop ((result != 0) ? 1 : 0);

			if(result)
			{
				nSharedAccessCount = 0;
			}

		}
	}

	if(result)
	{
		nExclusiveAccessCount++;
	}
	return result;
}

bool RWMutex::unlock()
{
	bool result1 = true, result2 = true;
	if(nExclusiveAccessCount == 0)
	{
		if(mtx_lock(&mtxSharedAccessCompleted) != thrd_success)
			return false;

//		nCompletedSharedAccessCount++;
		if(++nCompletedSharedAccessCount == 0)
		{
			result1 = (cnd_signal(&cndSharedAccessCompleted) == thrd_success);
		}
		result2 = (mtx_unlock(&mtxSharedAccessCompleted) == thrd_success);
	}
	else
	{
		nExclusiveAccessCount--;
		result1 = (mtx_unlock(&mtxSharedAccessCompleted) == thrd_success);
		result2 = (mtx_unlock(&mtxExclusiveAccess) == thrd_success);
	}

	return (result1 && result2);
}
