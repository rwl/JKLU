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

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;

/**
 * Sorts the columns of L and U so that the row indices appear in strictly
 * increasing order.
 *
 */
public class Dklu_sort extends Dklu_internal {

	/**
	 * Sort L or U using a double-transpose.
	 */
	public static void sort(int n, int[] Xip, int[] Xlen, double LU, int[] Tp,
			int[] Tj, double Tx, int[] W)
	{
		int[] Xi ;
		double Xx ;
		int p, i, j, len, nz, tp, xlen, pend ;

		ASSERT (Dklu_dump.klu_valid_LU (n, FALSE, Xip, Xlen, LU)) ;

		/* count the number of entries in each row of L or U */
		for (i = 0 ; i < n ; i++)
		{
			W [i] = 0 ;
		}
		for (j = 0 ; j < n ; j++)
		{
			GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
			for (p = 0 ; p < len ; p++)
			{
				W [Xi [p]]++ ;
			}
		}

		/* construct the row pointers for T */
		nz = 0 ;
		for (i = 0 ; i < n ; i++)
		{
			Tp [i] = nz ;
			nz += W [i] ;
		}
		Tp [n] = nz ;
		for (i = 0 ; i < n ; i++)
		{
			W [i] = Tp [i] ;
		}

		/* transpose the matrix into Tp, Ti, Tx */
		for (j = 0 ; j < n ; j++)
		{
			GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
			for (p = 0 ; p < len ; p++)
			{
				tp = W [Xi [p]]++ ;
				Tj [tp] = j ;
				Tx [tp] = Xx [p] ;
			}
		}

		/* transpose the matrix back into Xip, Xlen, Xi, Xx */
		for (j = 0 ; j < n ; j++)
		{
			W [j] = 0 ;
		}
		for (i = 0 ; i < n ; i++)
		{
			pend = Tp [i+1] ;
			for (p = Tp [i] ; p < pend ; p++)
			{
				j = Tj [p] ;
				GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
				xlen = W [j]++ ;
				Xi [xlen] = i ;
				Xx [xlen] = Tx [p] ;
			}
		}

		ASSERT (Dklu_dump.klu_valid_LU (n, FALSE, Xip, Xlen, LU)) ;
	}


	public static int klu_sort(KLU_symbolic Symbolic, KLU_numeric Numeric,
			KLU_common Common)
	{
		int[] R, W, Tp, Ti, Lip, Uip, Llen, Ulen ;
		double Tx ;
		double[][] LUbx ;
		int n, nk, nz, block, nblocks, maxblock, k1 ;
		int m1 ;

		if (Common == null)
		{
			return (FALSE) ;
		}
		Common.status = KLU_OK ;

		n = Symbolic.n ;
		R = Symbolic.R ;
		nblocks = Symbolic.nblocks ;
		maxblock = Symbolic.maxblock ;

		Lip  = Numeric.Lip ;
		Llen = Numeric.Llen ;
		Uip  = Numeric.Uip ;
		Ulen = Numeric.Ulen ;
		LUbx = (double[][]) Numeric.LUbx ;

		m1 = ((int) maxblock) + 1 ;

		/* allocate workspace */
		nz = MAX (Numeric.max_lnz_block, Numeric.max_unz_block) ;
		W  = KLU_malloc (maxblock, sizeof (Int), Common) ;
		Tp = KLU_malloc (m1, sizeof (Int), Common) ;
		Ti = KLU_malloc (nz, sizeof (Int), Common) ;
		Tx = KLU_malloc (nz, sizeof (double), Common) ;

		PRINTF ("\n======================= Start sort:\n") ;

		if (Common.status == KLU_OK)
		{
			/* sort each block of L and U */
			for (block = 0 ; block < nblocks ; block++)
			{
				k1 = R [block] ;
				nk = R [block+1] - k1 ;
				if (nk > 1)
				{
					PRINTF ("\n-------------------block: %d nk %d\n", block, nk) ;
					sort (nk, Lip + k1, Llen + k1, LUbx [block], Tp, Ti, Tx, W) ;
					sort (nk, Uip + k1, Ulen + k1, LUbx [block], Tp, Ti, Tx, W) ;
				}
			}
		}

		PRINTF ("\n======================= sort done.\n") ;

		/* free workspace */
		KLU_free (W, maxblock, sizeof (Int), Common) ;
		KLU_free (Tp, m1, sizeof (Int), Common) ;
		KLU_free (Ti, nz, sizeof (Int), Common) ;
		KLU_free (Tx, nz, sizeof (double), Common) ;
		return (Common.status == KLU_OK) ;
	}

}
