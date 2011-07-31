/**
 * KLU: a sparse LU factorization algorithm.
 * Copyright (C) 2004-2009, Timothy A. Davis.
 * Copyright (C) 2011, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/klu
 *
 * -------------------------------------------------------------------------
 *
 * KLU is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * KLU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

package edu.ufl.cise.klu.tdouble;

/**
 * KLU memory management routines.
 */
public class Dklu_memory {

	/**
	 * Safely compute a+b, and check for size_t overflow.
	 */
	public static size_t klu_add_size_t(size_t a, size_t b, int ok)
	{
	    ok = ok && ((a + b) >= max(a, b));
	    return ok ? (a + b) : (size_t) -1;
	}

	public static size_t klu_mult_size_t(size_t a, size_t k, Int ok)
	{
	    size_t i, s = 0;
	    for (i = 0 ; i < k ; i++)
	    {
	        s = klu_add_size_t(s, a, ok);
	    }
	    return ((ok) ? s : ((size_t) -1));
	}

	/**
	 * Wrapper around malloc routine (mxMalloc for a mexFunction).  Allocates
	 * space of size MAX(1,n)*size, where size is normally a sizeof (...).
	 *
	 * This routine and KLU_realloc do not set Common.status to KLU_OK on success,
	 * so that a sequence of KLU_malloc's or KLU_realloc's can be used.  If any of
	 * them fails, the Common.status will hold the most recent error status.
	 *
	 * Usage, for a pointer to Int:
	 *
	 *      p = KLU_malloc (n, sizeof (Int), Common)
	 *
	 * Uses a pointer to the malloc routine (or its equivalent) defined in Common.
	 *
	 * @param n number of items
	 * @param size size of each item
	 * @param Common
	 * @return
	 */
	public static Object klu_malloc(size_t n, size_t size, KLU_common Common)
	{
		Object p ;
	    size_t s ;
	    Int ok = true ;

	    if (Common == null)
	    {
	        p = null ;
	    }
	    else if (size == 0)
	    {
	        /* size must be > 0 */
	        Common.status = KLU_INVALID ;
	        p = null ;
	    }
	    else if (n >= INT_MAX)
	    {
	        /* object is too big to allocate; p[i] where i is an Int will not
	         * be enough. */
	        Common.status = KLU_TOO_LARGE ;
	        p = null ;
	    }
	    else
	    {
	        /* call malloc, or its equivalent */
	        s = klu_mult_size_t(max(1, n), size, ok) ;
	        p = ok ? ((Common.malloc_memory) (s)) : null ;
	        if (p == null)
	        {
	            /* failure: out of memory */
	            Common.status = KLU_OUT_OF_MEMORY ;
	        }
	        else
	        {
	            Common.memusage += s ;
	            Common.mempeak = max(Common.mempeak, Common.memusage) ;
	        }
	    }
	    return (p) ;
	}

	/**
	 * Wrapper around free routine.
	 *
	 * @param p block of memory to free
	 * @param n size of block to free, in # of items
	 * @param size size of each item
	 * @param Common
	 */
	public static void klu_free(Object p, size_t n, size_t size,
			KLU_common Common) {
		size_t s ;
	    Int ok = true ;
	    if (p != null && Common != null)
	    {
	        /* only free the object if the pointer is not null */
	        /* call free, or its equivalent */
	        (Common.free_memory) (p) ;
	        s = klu_mult_size_t(max(1, n), size, ok) ;
	        Common.memusage -= s ;
	    }
	    /* return null, and the caller should assign this to p.  This avoids
	     * freeing the same pointer twice. */
	    return (null) ;
	}

	/**
	 * Wrapper around realloc routine (mxRealloc for a mexFunction).  Given a
	 * pointer p to a block allocated by KLU_malloc, it changes the size of the
	 * block pointed to by p to be MAX(1,nnew)*size in size.  It may return a
	 * pointer different than p.  This should be used as (for a pointer to Int):
	 *
	 *      p = KLU_realloc (nnew, nold, sizeof (Int), p, Common) ;
	 *
	 * If p is null, this is the same as p = KLU_malloc (...).
	 * A size of nnew=0 is treated as nnew=1.
	 *
	 * If the realloc fails, p is returned unchanged and Common.status is set
	 * to KLU_OUT_OF_MEMORY.  If successful, Common.status is not modified,
	 * and p is returned (possibly changed) and pointing to a large block of memory.
	 *
	 * Uses a pointer to the realloc routine (or its equivalent) defined in Common.
	 *
	 * @param nnew requested # of items in reallocated block
	 * @param nold old # of items
	 * @param size size of each item
	 * @param p block of memory to realloc
	 * @param Common
	 * @return pointer to reallocated block
	 */
	public static Object klu_realloc(size_t nnew, size_t nold, size_t size,
			Object p, KLU_common Common)
	{
		Object pnew ;
	    size_t snew, sold ;
	    Int ok = true ;

	    if (Common == null)
	    {
	        p = null ;
	    }
	    else if (size == 0)
	    {
	        /* size must be > 0 */
	        Common.status = KLU_INVALID ;
	        p = null ;
	    }
	    else if (p == null)
	    {
	        /* A fresh object is being allocated. */
	        p = KLU_malloc (nnew, size, Common) ;
	    }
	    else if (nnew >= INT_MAX)
	    {
	        /* failure: nnew is too big.  Do not change p */
	        Common.status = KLU_TOO_LARGE ;
	    }
	    else
	    {
	        /* The object exists, and is changing to some other nonzero size. */
	        /* call realloc, or its equivalent */
	        snew = klu_mult_size_t(max(1,nnew), size, ok) ;
	        sold = klu_mult_size_t(max(1,nold), size, ok) ;
	        pnew = ok ? ((Common.realloc_memory) (p, snew)) : null ;
	        if (pnew == null)
	        {
	            /* Do not change p, since it still points to allocated memory */
	            Common.status = KLU_OUT_OF_MEMORY ;
	        }
	        else
	        {
	            /* success: return the new p and change the size of the block */
	            Common.memusage += (snew - sold) ;
	            Common.mempeak = max(Common.mempeak, Common.memusage) ;
	            p = pnew ;
	        }
	    }
	    return (p) ;
	}

}
