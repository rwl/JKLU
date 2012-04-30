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

import edu.ufl.cise.klu.common.KLU_common;

/**
 * Scale a matrix and check to see if it is valid.  Can be called by the user.
 * This is called by KLU_factor and KLU_refactor.  Returns true if the input
 * matrix is valid, false otherwise.  If the W input argument is non-null,
 * then the input matrix is checked for duplicate entries.
 *
 * scaling methods:
 *      <0: no scaling, do not compute Rs, and do not check input matrix.
 *      0: no scaling
 *      1: the scale factor for row i is sum (abs (A (i,:)))
 *      2 or more: the scale factor for row i is max(abs (A (i,:)))
 */
public class Dklu_scale extends Dklu_internal {

	/**
	 *
	 * @param scale 0: none, 1: sum, 2: max
	 * @param n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Ax outputs, not defined on input
	 * @param Rs size n, can be null if scale <= 0
	 * @param W size n, can be null
	 * @param Common
	 * @return true if successful, false otherwise
	 */
	public static int klu_scale(int scale, int n, int[] Ap, int[] Ai,
			double[] Ax, double[] Rs, int[] W, KLU_common Common)
	{
		double a ;
		double[] Az ;
		int row, col, p, pend ;
		boolean check_duplicates ;

		/* ---------------------------------------------------------------------- */
		/* check inputs */
		/* ---------------------------------------------------------------------- */

		if (Common == null)
		{
			return (FALSE) ;
		}
		Common.status = KLU_OK ;

		if (scale < 0)
		{
			/* return without checking anything and without computing the
			 * scale factors */
			return (TRUE) ;
		}

		Az = (double[]) Ax ;

		if (n <= 0 || Ap == null || Ai == null || Az == null ||
			(scale > 0 && Rs == null))
		{
			/* Ap, Ai, Ax and Rs must be present, and n must be > 0 */
			Common.status = KLU_INVALID ;
			return (FALSE) ;
		}
		if (Ap [0] != 0 || Ap [n] < 0)
		{
			/* nz = Ap [n] must be >= 0 and Ap [0] must equal zero */
			Common.status = KLU_INVALID ;
			return (FALSE) ;
		}
		for (col = 0 ; col < n ; col++)
		{
			if (Ap [col] > Ap [col+1])
			{
				/* column pointers must be non-decreasing */
				Common.status = KLU_INVALID ;
				return (FALSE) ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* scale */
		/* ---------------------------------------------------------------------- */

		if (scale > 0)
		{
			/* initialize row sum or row max */
			for (row = 0 ; row < n ; row++)
			{
				Rs [row] = 0 ;
			}
		}

		/* check for duplicates only if W is present */
		check_duplicates = (W != null) ;
		if (check_duplicates)
		{
			for (row = 0 ; row < n ; row++)
			{
				W [row] = EMPTY ;
			}
		}

		for (col = 0 ; col < n ; col++)
		{
			pend = Ap [col+1] ;
			for (p = Ap [col] ; p < pend ; p++)
			{
				row = Ai [p] ;
				if (row < 0 || row >= n)
				{
					/* row index out of range, or duplicate entry */
					Common.status = KLU_INVALID ;
					return (FALSE) ;
				}
				if (check_duplicates)
				{
					if (W [row] == col)
					{
						/* duplicate entry */
						Common.status = KLU_INVALID ;
						return (FALSE) ;
					}
					/* flag row i as appearing in column col */
					W [row] = col ;
				}
				a = ABS (Az [p]) ;
				//ABS (a, Az [p]) ;
				if (scale == 1)
				{
					/* accumulate the abs. row sum */
					Rs [row] += a ;
				}
				else if (scale > 1)
				{
					/* find the max abs. value in the row */
					Rs [row] = MAX (Rs [row], a) ;
				}
			}
		}

		if (scale > 0)
		{
			/* do not scale empty rows */
			for (row = 0 ; row < n ; row++)
			{
				/* matrix is singular */
				PRINTF ("Rs [%d] = %g\n", row, Rs [row]) ;

				if (Rs [row] == 0.0)
				{
					PRINTF ("Row %d of A is all zero\n", row) ;
					Rs [row] = 1.0 ;
				}
			}
		}

		return (TRUE) ;
	}

}
