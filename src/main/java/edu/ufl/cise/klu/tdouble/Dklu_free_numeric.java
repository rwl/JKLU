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

/**
 * Free the KLU Numeric object.
 */
public class Dklu_free_numeric extends Dklu_internal
{

//	public static int klu_free_numeric(KLU_numeric NumericHandle,
//			KLU_common  Common)
//	{
//		KLU_numeric Numeric ;
//		double[][] LUbx ;
//		int LUsize ;
//		int block, n, nzoff, nblocks ;
//
//		if (Common == null)
//		{
//			return (FALSE) ;
//		}
//		if (NumericHandle == null || NumericHandle == null)
//		{
//			return (TRUE) ;
//		}
//
//		Numeric = NumericHandle ;
//
//		n = Numeric.n ;
//		nzoff = Numeric.nzoff ;
//		nblocks = Numeric.nblocks ;
//		LUsize = Numeric.LUsize ;
//
//		LUbx = (double[][]) Numeric.LUbx ;
//		if (LUbx != null)
//		{
//			for (block = 0 ; block < nblocks ; block++)
//			{
//				KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
//					sizeof (double), Common) ;
//			}
//		}
//
//		KLU_free (Numeric.Pnum, n, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Offp, n+1, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Offi, nzoff+1, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Offx, nzoff+1, sizeof (double), Common) ;
//
//		KLU_free (Numeric.Lip,  n, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Llen, n, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Uip,  n, sizeof (Integer), Common) ;
//		KLU_free (Numeric.Ulen, n, sizeof (Integer), Common) ;
//
//		KLU_free (Numeric.LUsize, nblocks, sizeof (int), Common) ;
//
//		KLU_free (Numeric.LUbx, nblocks, sizeof (double[]), Common) ;
//
//		KLU_free (Numeric.Udiag, n, sizeof (double), Common) ;
//
//		KLU_free (Numeric.Rs,   n, sizeof (Double), Common) ;
//		KLU_free (Numeric.Pinv, n, sizeof (Int), Common) ;
//
//		KLU_free (Numeric.Work, Numeric.worksize, 1, Common) ;
//
//		KLU_free (Numeric, 1, sizeof (KLU_numeric), Common) ;
//
//		NumericHandle = null ;
//		return (TRUE) ;
//	}

}
