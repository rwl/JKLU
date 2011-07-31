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

public class Dklu_version {

	/** enable diagnostic printing */
	static boolean NPRINT = true;

	static final int INT_MAX = 0x7fffffff;

	static final String INT_ID = "%d";

	static int BYTES (Object type, double n)
	{
		return sizeof (type * n);
	}

	static double CEILING (double b, double u)
	{
		return (b+u-1) / u;
	}

	static double UNITS (Object type, double n)
	{
		return CEILING (BYTES (type, n), sizeof (Unit));
	}

	static double DUNITS (Object type, int n)
	{
		return Math.ceil(BYTES (type, (double) n) / sizeof (Unit));
	}

	static void GET_I_POINTER(LU, Xip, Xi, k)
	{
	    Xi = (Int[]) (LU + Xip [k]) ;
	}

	static void GET_X_POINTER(LU, Xip, Xlen, Xx, k)
	{
	    Xx = (Entry[]) (LU + Xip [k] + UNITS (Int, Xlen [k])) ;
	}

	static void GET_POINTER(LU, Xip, Xlen, Xi, Xx, k, xlen)
	{
	    Unit xp = LU + Xip [k] ;
	    xlen = Xlen [k] ;
	    Xi = (Int[]) xp ;
	    Xx = (Entry[]) (xp + UNITS (Int, xlen)) ;
	}

	static boolean SCALAR_IS_NAN (double x)
	{
		return x != x;
	}

	static boolean SCALAR_IS_ZERO (double x)
	{
		return x == 0.0;
	}

	static boolean SCALAR_IS_NONZERO (double x)
	{
		return x != 0.0;
	}

	static boolean SCALAR_IS_LTZERO (double x)
	{
		return x < 0.0;
	}

	/* scalar absolute value macro. If x is NaN, the result is NaN */
	static double SCALAR_ABS (double x)
	{
		return SCALAR_IS_LTZERO (x) ? -x : x;
	}

	static void PRINTF (String format, Object... args)
	{
		System.out.printf(format, args);
	}

	static double PRINT_SCALAR (double a)
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

	static final Class Unit = Double.class;

	static final Class Entry = Double.class;

	/**
	 * @return TRUE if a complex number is in split form, FALSE if in packed
	 * form.
	 */
	static int SPLIT (double s)
	{
		return 1 ;
	}

	/**
	 * @return real part of c
	 */
	static double REAL (double c)
	{
		return c ;
	}

	/**
	 * @return imag part of c
	 */
	static double IMAG (double c)
	{
		return 0.0 ;
	}

	/**
	 * c = (s1) + (s2)*i
	 */
	static void ASSIGN (Double c, double[] s1, double[] s2, int p,
			boolean split)
	{
		c = s1[p] ;
	}

	static void CLEAR (Double c)
	{
		c = 0.0 ;
	}

	static void CLEAR_AND_INCREMENT (Double p)
	{
		p = 0.0 ;
		p++;
	}

	/**
	 * @return True if a is NaN
	 */
	static boolean IS_NAN (double a)
	{
		return SCALAR_IS_NAN (a) ;
	}

	/**
	 * @return True if a == 0
	 */
	static boolean IS_ZERO (double a)
	{
		return SCALAR_IS_ZERO (a) ;
	}

	/**
	 * @return True if a != 0
	 */
	static boolean IS_NONZERO (double a)
	{
		return SCALAR_IS_NONZERO (a) ;
	}

	/**
	 * c /= s
	 */
	static void SCALE_DIV (Double c, double s)
	{
		c /= s ;
	}

	/**
	 * a = c/s
	 */
	static void SCALE_DIV_ASSIGN (Double a, double c, double s)
	{
		a = c / s ;
	}

	/**
	 * c *= s
	 */
	static void SCALE (Double c, double s)
	{
		c *= s ;
	}

	/**
	 * c += a
	 */
	static void ASSEMBLE (Double c, double a)
	{
		c += a ;
	}

	/**
	 * c += *p++
	 */
	static void ASSEMBLE_AND_INCREMENT (Double c, double p)
	{
		c += p++ ;
	}

	/**
	 * c -= a
	 */
	static void DECREMENT (Double c, double a)
	{
		c -= a ;
	}

	/**
	 * c = a*b
	 */
	static void MULT (Double c, double a,  double b)
	{
		c = a * b ;
	}

	/**
	 * c = a*conjugate(b)
	 */
	static void MULT_CONJ (Double c, double a, double b)
	{
		c = a * b ;
	}

	/**
	 * c -= a*b
	 */
	static void MULT_SUB (Double c, double a, double b)
	{
		c -= a * b ;
	}

	/**
	 * c -= a*conjugate(b)
	 */
	static void MULT_SUB_CONJ (Double c, double a, double b)
	{
		c -= a * b ;
	}

	/**
	 * c = a/b
	 */
	static void DIV (Double c, double a, double b)
	{
		c = a / b ;
	}

	/**
	 * c = 1/c
	 */
	static void RECIPROCAL (Double c)
	{
		c = 1.0 / c ;
	}

	/**
	 * c = a/conjugate(b)
	 */
	static void DIV_CONJ (Double c, double a, double b)
	{
		c = a / b ;
	}

	/**
	 * approximate absolute value, s = |r|+|i|
	 */
	static void APPROX_ABS (Double s, double a)
	{
		s = SCALAR_ABS (a) ;
	}

	/**
	 * exact absolute value, s = sqrt (a.real^2 + amag^2)
	 */
	static void ABS (Double s, double a)
	{
		s = SCALAR_ABS (a) ;
	}

	static void PRINT_ENTRY (double a)
	{
		PRINT_SCALAR (a) ;
	}

	static void CONJ (Double a, double x)
	{
		a = x;
	}

	/* for flop counts */
	static final double MULTSUB_FLOPS   = 2.0 ;      /* c -= a*b */
	static final double DIV_FLOPS       = 1.0 ;      /* c = a/b */
	static final double ABS_FLOPS       = 0.0 ;      /* c = abs (a) */
	static final double ASSEMBLE_FLOPS  = 1.0 ;      /* c += a */
	static final double DECREMENT_FLOPS = 1.0 ;      /* c -= a */
	static final double MULT_FLOPS      = 1.0 ;      /* c = a*b */
	static final double SCALE_FLOPS     = 1.0 ;      /* c = a/s */

}
