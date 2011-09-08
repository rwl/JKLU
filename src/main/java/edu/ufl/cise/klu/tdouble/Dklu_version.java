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

public abstract class Dklu_version {

	public static final int KLU_OK = 0;
	/** status > 0 is a warning, not an error */
	public static final int KLU_SINGULAR = 1;
	public static final int KLU_OUT_OF_MEMORY = -2;
	public static final int KLU_INVALID = -3;
	/** integer overflow has occured */
	public static final int KLU_TOO_LARGE = -4;

	/** enable diagnostic printing */
	protected static boolean NPRINT = true ;

	protected static final int INT_MAX = 0x7fffffff ;

	protected static final String INT_ID = "%d" ;

	protected static int BYTES (Object type, double n)
	{
		return sizeof (type * n) ;
	}

	protected static double CEILING (double b, double u)
	{
		return (b+u-1) / u ;
	}

	protected static double UNITS (Object type, double n)
	{
		return CEILING (BYTES (type, n), sizeof (double)) ;
	}

	protected static double DUNITS (Object type, int n)
	{
		return Math.ceil(BYTES (type, (double) n) / sizeof (double)) ;
	}

	protected static void GET_I_POINTER(LU, Xip, Xi, k)
	{
		Xi = (Int[]) (LU + Xip [k]) ;
	}

	protected static void GET_X_POINTER(double[] LU, int[] Xip, int Xlen,
			double[] Xx, int k)
	{
		Xx = (double[]) (LU + Xip [k] + UNITS (Int, Xlen [k])) ;
	}

	protected static void GET_POINTER(double[] LU, int[] Xip, int Xlen,
			int[] Xi, double[] Xx, int k, int xlen)
	{
		double xp = LU + Xip [k] ;
		xlen = Xlen [k] ;
		Xi = (Int[]) xp ;
		Xx = (double[]) (xp + UNITS (Int, xlen)) ;
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
		System.out.printf(format, args) ;
	}

	protected static double PRINT_SCALAR (double a)
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
	/* Real floating-point arithmetic */
	/* ---------------------------------------------------------------------- */

	/**
	 * @return TRUE if a complex number is in split form, FALSE if in packed
	 * form.
	 */
	protected static int SPLIT (double s)
	{
		return 1 ;
	}

	/**
	 * @return real part of c
	 */
	protected static double REAL (double c)
	{
		return c ;
	}

	/**
	 * @return imag part of c
	 */
	protected static double IMAG (double c)
	{
		return 0.0 ;
	}

	/**
	 * c = (s1) + (s2)*i
	 */
//	protected static void ASSIGN (Double c, double[] s1, double[] s2, int p,
//			boolean split)
//	{
//		c = s1[p] ;
//	}

//	protected static void CLEAR (Double c)
//	{
//		c = 0.0 ;
//	}

//	protected static void CLEAR_AND_INCREMENT (Double p)
//	{
//		p = 0.0 ;
//		p++ ;
//	}

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

	/**
	 * @return True if a != 0
	 */
	protected static boolean IS_NONZERO (double a)
	{
		return SCALAR_IS_NONZERO (a) ;
	}

	/**
	 * c /= s
	 */
//	protected static void SCALE_DIV (Double c, double s)
//	{
//		c /= s ;
//	}

	/**
	 * a = c/s
	 */
	protected static double SCALE_DIV_ASSIGN (double c, double s)
	{
		return c / s ;
	}

	/**
	 * c *= s
	 */
//	protected static void SCALE (Double c, double s)
//	{
//		c *= s ;
//	}

	/**
	 * c += a
	 */
//	protected static void ASSEMBLE (Double c, double a)
//	{
//		c += a ;
//	}

	/**
	 * c += *p++
	 */
//	protected static void ASSEMBLE_AND_INCREMENT (Double c, double p)
//	{
//		c += p++ ;
//	}

	/**
	 * c -= a
	 */
//	protected static void DECREMENT (Double c, double a)
//	{
//		c -= a ;
//	}

	/**
	 * c = a*b
	 */
//	protected static void MULT (Double c, double a,  double b)
//	{
//		c = a * b ;
//	}

	/**
	 * c = a*conjugate(b)
	 */
//	protected static void MULT_CONJ (Double c, double a, double b)
//	{
//		c = a * b ;
//	}

	/**
	 * c -= a*b
	 */
//	protected static void MULT_SUB (Double c, double a, double b)
//	{
//		c -= a * b ;
//	}
//
//	/**
//	 * c -= a*conjugate(b)
//	 */
//	protected static void MULT_SUB_CONJ (Double c, double a, double b)
//	{
//		c -= a * b ;
//	}
//
//	/**
//	 * c = a/b
//	 */
//	protected static void DIV (Double c, double a, double b)
//	{
//		c = a / b ;
//	}
//
//	/**
//	 * c = 1/c
//	 */
//	protected static void RECIPROCAL (Double c)
//	{
//		c = 1.0 / c ;
//	}
//
//	/**
//	 * c = a/conjugate(b)
//	 */
//	protected static void DIV_CONJ (Double c, double a, double b)
//	{
//		c = a / b ;
//	}
//
//	/**
//	 * approximate absolute value, s = |r|+|i|
//	 */
//	protected static void APPROX_ABS (Double s, double a)
//	{
//		s = SCALAR_ABS (a) ;
//	}
//
//	/**
//	 * exact absolute value, s = sqrt (a.real^2 + amag^2)
//	 */
//	protected static void ABS (Double s, double a)
//	{
//		s = SCALAR_ABS (a) ;
//	}

	protected static void PRINT_ENTRY (double a)
	{
		PRINT_SCALAR (a) ;
	}

//	protected static void CONJ (Double a, double x)
//	{
//		a = x ;
//	}

	/* for flop counts */
	protected static final double MULTSUB_FLOPS   = 2.0 ;      /* c -= a*b */
	protected static final double DIV_FLOPS       = 1.0 ;      /* c = a/b */
	protected static final double ABS_FLOPS       = 0.0 ;      /* c = abs(a) */
	protected static final double ASSEMBLE_FLOPS  = 1.0 ;      /* c += a */
	protected static final double DECREMENT_FLOPS = 1.0 ;      /* c -= a */
	protected static final double MULT_FLOPS      = 1.0 ;      /* c = a*b */
	protected static final double SCALE_FLOPS     = 1.0 ;      /* c = a/s */

}
