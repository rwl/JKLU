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

package edu.ufl.cise.klu.tdouble;

public abstract class Dklu_internal extends Dklu_version {

	/**
	 * enable debugging and assertions
	 */
	public static boolean NDEBUG = true ;

	protected static void ASSERT (boolean a)
	{
		if (!NDEBUG)
		{
			assert a ;
		}
	}

	protected static void ASSERT (int a)
	{
		ASSERT (a != 0) ;
	}

	/**
	 * @return true if an integer (stored in double x) would overflow (or if
	 * x is NaN)
	 */
	protected static boolean INT_OVERFLOW (double x)
	{
		return ((!(x * (1.0+1e-8) <= (double) INT_MAX))
							|| SCALAR_IS_NAN (x)) ;
	}

	protected static final int TRUE = 1 ;
	protected static final int FALSE = 0 ;

	protected static int MAX (int a, int b)
	{
		return a > b ?  a : b ;
	}

	protected static int MIN (int a, int b)
	{
		return a < b ?  a : b ;
	}

	protected static double MAX (double a, double b)
	{
		return a > b ?  a : b ;
	}

	protected static double MIN (double a, double b)
	{
		return a < b ?  a : b ;
	}

	protected static long MAX (long a, long b)
	{
		return a > b ?  a : b ;
	}

	/* FLIP is a "negation about -1", and is used to mark an integer i that is
	 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
	 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
	 * for all integers i.  UNFLIP (i) is >= EMPTY. */
	protected static final int EMPTY = -1 ;

	protected static int FLIP (int i)
	{
		return -i - 2 ;
	}

	protected static double FLIP (double i)
	{
		return -i - 2 ;
	}

	protected static int UNFLIP (int i)
	{
		return (i < EMPTY) ? FLIP (i) : i ;
	}

	protected static double UNFLIP (double i)
	{
		return (i < EMPTY) ? FLIP (i) : i ;
	}

}
