/**
 * KLU: a sparse LU factorization algorithm.
 * Copyright (C) 2004-2009, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
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

package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;

/**
 * KLU memory management routines.
 */
public class DZklu_memory extends DZklu_internal {

	/**
	 * Safely compute a+b, and check for int overflow.
	 */
	public static int klu_z_add_size_t(int a, int b, int[] ok)
	{
		(ok[0]) = (ok[0] != 0) && ((a + b) >= MAX (a,b)) ? 1 : 0;
		return ((ok[0] != 0) ? (a + b) : ((int) -1)) ;
	}

	public static int klu_z_mult_size_t(int a, int k, int[] ok)
	{
		int i, s = 0 ;
		for (i = 0 ; i < k ; i++)
		{
			s = klu_z_add_size_t (s, a, ok) ;
		}
		return ((ok[0] != 0) ? s : ((int) -1)) ;
	}

	/**
	 * Allocates space of size MAX(1,n).
	 *
	 * This routine and KLU_realloc do not set Common.status to KLU_OK on success,
	 * so that a sequence of KLU_malloc's or KLU_realloc's can be used.  If any of
	 * them fails, the Common.status will hold the most recent error status.
	 *
	 * Usage, for a pointer to Int:
	 *
	 *      p = KLU_malloc (n, sizeof (Int), Common)
	 *
	 * @param n number of items
	 * @param size size of each item
	 * @param Common
	 * @return
	 */
	public static int[] klu_z_malloc_int(int n, KLU_common Common)
	{
		Runtime runtime;
		int[] p = null;

		if (n >= INT_MAX)
		{
			Common.status = KLU_TOO_LARGE ;
			p = null ;
		}
		else
		{
			try
			{
				p = new int[n];
				runtime = Runtime.getRuntime ();
				Common.memusage = runtime.totalMemory () - runtime.freeMemory ();
				Common.mempeak = MAX (Common.mempeak, Common.memusage) ;
			}
			catch (OutOfMemoryError e)
			{
				/* failure: out of memory */
				Common.status = KLU_OUT_OF_MEMORY ;
				p = null;
			}
		}
		return (p) ;
	}

	public static double[] klu_z_malloc_dbl(int n, KLU_common Common)
	{
		Runtime runtime;
		double[] p = null;

		if (n >= INT_MAX)
		{
			Common.status = KLU_TOO_LARGE ;
			p = null ;
		}
		else
		{
			try
			{
				p = new double[n];
				runtime = Runtime.getRuntime ();
				Common.memusage = runtime.totalMemory () - runtime.freeMemory ();
				Common.mempeak = MAX (Common.mempeak, Common.memusage) ;
			}
			catch (OutOfMemoryError e)
			{
				/* failure: out of memory */
				Common.status = KLU_OUT_OF_MEMORY ;
				p = null;
			}
		}
		return (p) ;
	}

	/**
	 * Given an array p allocated by KLU_malloc, it changes the size of the
	 * block pointed to by p to be MAX(1,nnew) in size.  It may return an
	 * array different than p.  This should be used as:
	 *
	 *      p = KLU_realloc (nnew, nold, p, Common) ;
	 *
	 * If p is null, this is the same as p = KLU_malloc (...).
	 * A size of nnew=0 is treated as nnew=1.
	 *
	 * If the realloc fails, p is returned unchanged and Common.status is set
	 * to KLU_OUT_OF_MEMORY.  If successful, Common.status is not modified,
	 * and p is returned (possibly changed) and pointing to a large block of memory.
	 *
	 * @param nnew requested # of items in reallocated block
	 * @param nold old # of items
	 * @param p block of memory to realloc
	 * @param Common
	 * @return pointer to reallocated block
	 */
	public static double[] klu_z_realloc_dbl (int nnew, int nold,
			double[] p, KLU_common Common)
	{
		double[] pnew ;
		int snew ;
		int sold ;

		if (Common == null)
		{
			p = null ;
		}
		else if (p == null)
		{
			/* A fresh object is being allocated. */
			p = klu_z_malloc_dbl (nnew, Common) ;
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
			snew = MAX (1, nnew) ;
			sold = MAX (1, nold) ;
			try
			{
				pnew = new double[snew] ;
				System.arraycopy(p, 0, pnew, 0, MIN (snew, sold)) ;
				Runtime runtime = Runtime.getRuntime();
				Common.memusage = runtime.totalMemory() - runtime.freeMemory();
				//Common.memusage += (snew - sold) ;
				Common.mempeak = MAX (Common.mempeak, Common.memusage) ;
				p = pnew ;
			}
			catch (OutOfMemoryError e)
			{
				/* Do not change p, since it still points to allocated memory */
				Common.status = KLU_OUT_OF_MEMORY ;
			}
		}
		return (p) ;
	}

}
