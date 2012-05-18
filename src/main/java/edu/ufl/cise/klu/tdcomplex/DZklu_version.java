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

import edu.ufl.cise.klu.tdcomplex.DZklu_common.DZklua;

public abstract class DZklu_version {

	public static final int KLU_OK = 0;
	/** status > 0 is a warning, not an error */
	public static final int KLU_SINGULAR = 1;
	public static final int KLU_OUT_OF_MEMORY = -2;
	public static final int KLU_INVALID = -3;
	/** integer overflow has occured */
	public static final int KLU_TOO_LARGE = -4;

	/** enable diagnostic printing */
	public static boolean NPRINT = true ;

	protected static final int INT_MAX = 0x7fffffff ;

	protected static final String INT_ID = "%d" ;

	public static final double[] CZERO = new double[] {0.0, 0.0} ;
	public static final double[] CONE = new double[] {1.0, 0.0} ;

	protected static double[] GET_I_POINTER(double[] LU, int[] Xip,
			int Xip_offset, int[] Xi_offset, int k)
	{
		Xi_offset[0] = Xip [Xip_offset + k] ;
		return LU ;
	}

	protected static double[] GET_POINTER(double[] LU,
			int[] Xip, int Xip_offset,
			int[] Xlen, int Xlen_offset,
			int[] Xi_offset,
			int[] Xx_offset,
			int k, int[] xlen)
	{
		int xp = Xip [Xip_offset + k] ;
		xlen[0] = Xlen [Xlen_offset + k] ;
		Xi_offset[0] = xp ;
		Xx_offset[0] = xp + xlen[0] ;
		return LU ;
	}

	protected static boolean SCALAR_IS_NAN (double x)
	{
		return x != x ;
	}

	protected static boolean SCALAR_IS_ZERO (double x)
	{
		return x == 0.0 ;
	}

	protected static boolean SCALAR_IS_NONZERO (double x)
	{
		return x != 0.0 ;
	}

	protected static boolean SCALAR_IS_LTZERO (double x)
	{
		return x < 0.0 ;
	}

	/* scalar absolute value macro. If x is NaN, the result is NaN */
	protected static double SCALAR_ABS (double x)
	{
		return SCALAR_IS_LTZERO (x) ? -x : x ;
	}

	protected static void PRINTF (String format, Object... args)
	{
		if (!NPRINT) {
			System.out.printf(format, args) ;
		}
	}

	protected static void PRINT_SCALAR (double a)
	{
		if (!NPRINT) {

			if (SCALAR_IS_NONZERO (a))
			{
				PRINTF (" (%g)", a) ;
			}
			else
			{
				PRINTF (" (0)") ;
			}

		}
	}

	/* ---------------------------------------------------------------------- */
	/* Complex floating-point arithmetic */
	/* ---------------------------------------------------------------------- */

	/**
	 * @return real part of c at index idx
	 */
	protected static double REAL (double[] c, int idx)
	{
		return c [2 * idx] ;
	}

	/**
	 * @return imag part of c at index idx
	 */
	protected static double IMAG (double[] c, int idx)
	{
		return c [(2 * idx) + 1] ;
	}

	/**
	 * @return True if a is NaN
	 */
	protected static boolean IS_NAN (double a)
	{
		return SCALAR_IS_NAN (a) ;
	}

	/**
	 * @return True if a == 0
	 */
	protected static boolean IS_ZERO (double a)
	{
		return SCALAR_IS_ZERO (a) ;
	}

	protected static boolean IS_ZERO (double[] a)
	{
		return SCALAR_IS_ZERO (a [0]) && SCALAR_IS_ZERO (a [1]) ;
	}

	public static final double [] CPLUS(double [] x, double [] y)
	{
		return new double [] {x [0] + y [0], x [1] + y [1]} ;
	}

	public static final double [] CPLUS(double [] x, double real)
	{
		return new double [] {x [0] + real, x [1]} ;
	}

	protected static boolean IS_NONZERO (double[] a)
	{
		return SCALAR_IS_NONZERO (a [0]) || SCALAR_IS_NONZERO (a [0]) ;
	}

	protected static void SCALE_DIV (double[] c, double s)
	{
		c [0] /= s ;
		c [1] /= s ;
	}

	/**
	 * a = c/s
	 */
	protected static void SCALE_DIV_ASSIGN (double[] a, double[] c, double s)
	{
		a [0] = c [0] / s ;
		a [1] = c [1] / s ;
	}

	protected static void SCALE_DIV_ASSIGN (DZklua c, int idx, double[] b, double s)
	{
		c.real(idx, b [0] / s) ;
		c.imag(idx, b [1] / s) ;
	}

	protected static void DIV (DZklua a, int idx, double[] b, double[] c) {
		a.set(idx, DIV(b, c));
	}

	/* This uses ACM Algo 116, by R. L. Smith, 1962. */
	/* c can be the same variable as a or b. */
	/* Ignore NaN case for double relop br>=bi. */
	protected static double[] DIV (double[] a, double[] b) {
		double r, den, ar, ai, br, bi ;
		double[] c ;
		br = b [0] ;
		bi = b [1] ;
		ar = a [0] ;
		ai = a [1] ;
		if (SCALAR_ABS (br) >= SCALAR_ABS (bi))
		{
			r = bi / br ;
			den = br + r * bi ;
			c = new double[] {
				(ar + ai * r) / den,
				(ai - ar * r) / den
			} ;
		}
		else
		{
			r = br / bi ;
			den = r * br + bi ;
			c = new double[] {
				(ar * r + ai) / den,
				(ai * r - ar) / den
			} ;
		}
		return c ;
	}

	protected static void MULT_SUB_CONJ (DZklua c, int idx, double[] a, double[] b)
	{
		double[] z = c.get(idx) ;
		MULT_SUB_CONJ (z, a, b) ;
		c.set(idx, z) ;
	}

	protected static void MULT_SUB_CONJ (double[] c, double[] a, double[] b)
	{
		c [0] -= a [0] * b [0] + a [1] * b [1] ;
		c [1] -= a [1] * b [0] - a [0] * b [1] ;
	}

	protected static void MULT_SUB (DZklua c, int idx, double[] a, double[] b)
	{
		double[] z = c.get(idx) ;
		MULT_SUB (z, a, b) ;
		c.set(idx, z) ;
	}

	protected static void MULT_SUB (double[] c, double[] a, double[] b)
	{
		c [0] -= a [0] * b [0] - a [1] * b [1] ;
		c [1] -= a [1] * b [0] + a [0] * b [1] ;
	}

	protected static void PRINT_ENTRY (double[] a)
	{
		if (SCALAR_IS_NONZERO (a [0]))
		{
			PRINTF (" (%g", a [0]) ;
		}
		else
		{
			PRINTF (" (0") ;
		}

		if (SCALAR_IS_LTZERO (a [1]))
		{
			PRINTF (" - %gi)", -a [1]) ;
		}
		else if (SCALAR_IS_ZERO (a [1]))
		{
			PRINTF (" + 0i)") ;
		}
		else
		{
			PRINTF (" + %gi)", a [1]) ;
		}
	}

	protected static double ABS (double[] a)
	{
		double r, ar, ai, s ;
		ar = SCALAR_ABS (a [0]) ;
		ai = SCALAR_ABS (a [1]) ;
		if (ar >= ai)
		{
			if (ar + ai == ar)
			{
				(s) = ar ;
			}
			else
		        {
				r = ai / ar ;
				(s) = ar * Math.sqrt (1.0 + r*r) ;
		        }
		}
		else
		{
			if (ai + ar == ai)
		        {
				(s) = ai ;
		        }
		        else
		        {
		        	r = ar / ai ;
		        	(s) = ai * Math.sqrt (1.0 + r*r) ;
		        }
		}
		return (s) ;
	}

	protected static void CLEAR(DZklua A, int i)
	{
		A.set(i, CZERO) ;
	}

	protected static double[] CONJ(double[] a)
	{
		return new double[] {a [0], -a [1]} ;
	}

	protected static void CLEAR(double[] A, int i)
	{
		A [i] = 0.0 ;
	}

	protected static void CLEAR(double[] a)
	{
		a [0] = 0.0 ;
		a [1] = 0.0 ;
	}

	/* for flop counts */
	protected static final double MULTSUB_FLOPS   = 8.0 ;      /* c -= a*b */
	protected static final double DIV_FLOPS       = 9.0 ;      /* c = a/b */
	protected static final double ABS_FLOPS       = 6.0 ;      /* c = abs(a) */
	protected static final double ASSEMBLE_FLOPS  = 2.0 ;      /* c += a */
	protected static final double DECREMENT_FLOPS = 2.0 ;      /* c -= a */
	protected static final double MULT_FLOPS      = 6.0 ;      /* c = a*b */
	protected static final double SCALE_FLOPS     = 2.0 ;      /* c = a/s */

}
