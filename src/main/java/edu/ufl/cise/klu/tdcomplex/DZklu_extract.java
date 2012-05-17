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
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;

/**
 * Extract KLU factorization into conventional compressed-column matrices.
 * If any output array is null, that part of the LU factorization is not
 * extracted (this is not an error condition).
 *
 * nnz(L) = Numeric.lnz, nnz(U) = Numeric.unz, and nnz(F) = Numeric.Offp [n]
 */
public class DZklu_extract extends DZklu_internal
{

	/**
	 *
	 * @param Numeric
	 * @param Symbolic
	 * @param Lp size n+1
	 * @param Li size nnz(L)
	 * @param Lx size nnz(L)
	 * @param Up size n+1
	 * @param Ui size nnz(U)
	 * @param Ux size nnz(U)
	 * @param Fp size n+1
	 * @param Fi size nnz(F)
	 * @param Fx size nnz(F)
	 * @param P row permutation, size n
	 * @param Q column permutation, size n
	 * @param Rs scale factors, size n
	 * @param R block boundaries, size nblocks+1
	 * @param Common
	 * @return
	 */
	public static int klu_z_extract(KLU_numeric Numeric, KLU_symbolic Symbolic,
			int[] Lp, int[] Li, double[] Lx, int[] Up, int[] Ui, double[] Ux,
			int[] Fp, int[] Fi, double[] Fx, int[] P, int[] Q, double[] Rs,
			int[] R, KLU_common Common)
	{
		int[] Lip, Llen, Uip, Ulen ;
		/*int[]*/double[] Li2, Ui2 ;
		double[] LU ;
		double[] Lx2, Ux2, Ukk ;
		int i, k, block, nblocks, n, nz, k1, k2, nk, kk, p ;
		int[] len = new int[1] ;
		int[] Li2_offset = new int[1] ;
		int[] Lx2_offset = new int[1] ;
		int[] Ui2_offset = new int[1] ;
		int[] Ux2_offset = new int[1] ;

		if (Common == null)
		{
			return (FALSE) ;
		}

		if (Symbolic == null || Numeric == null)
		{
			Common.status = KLU_INVALID ;
			return (FALSE) ;
		}

		Common.status = KLU_OK ;
		n = Symbolic.n ;
		nblocks = Symbolic.nblocks ;

		/* ---------------------------------------------------------------------- */
		/* extract scale factors */
		/* ---------------------------------------------------------------------- */

		if (Rs != null)
		{
			if (Numeric.Rs != null)
			{
				for (i = 0 ; i < n ; i++)
				{
					Rs [i] = Numeric.Rs [i] ;
				}
			}
			else
			{
				/* no scaling */
				for (i = 0 ; i < n ; i++)
				{
					Rs [i] = 1 ;
				}
			}
		}

		/* ---------------------------------------------------------------------- */
		/* extract block boundaries */
		/* ---------------------------------------------------------------------- */

		if (R != null)
		{
			for (block = 0 ; block <= nblocks ; block++)
			{
				R [block] = Symbolic.R [block] ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* extract final row permutation */
		/* ---------------------------------------------------------------------- */

		if (P != null)
		{
			for (k = 0 ; k < n ; k++)
			{
				P [k] = Numeric.Pnum [k] ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* extract column permutation */
		/* ---------------------------------------------------------------------- */

		if (Q != null)
		{
			for (k = 0 ; k < n ; k++)
			{
				Q [k] = Symbolic.Q [k] ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* extract each block of L */
		/* ---------------------------------------------------------------------- */

		if (Lp != null && Li != null && Lx != null)
		{
			nz = 0 ;
			for (block = 0 ; block < nblocks ; block++)
			{
				k1 = Symbolic.R [block] ;
				k2 = Symbolic.R [block+1] ;
				nk = k2 - k1 ;
				if (nk == 1)
				{
					/* singleton block */
					Lp [k1] = nz ;
					Li [nz] = k1 ;
					Lx [nz] = 1 ;
					nz++ ;
				}
				else
				{
					/* non-singleton block */
					LU = Numeric.LUbx [block] ;
					Lip = Numeric.Lip ;
					int Lip_offset = k1 ;
					Llen = Numeric.Llen ;
					int Llen_offset = k1 ;
					for (kk = 0 ; kk < nk ; kk++)
					{
						Lp [k1+kk] = nz ;
						/* add the unit diagonal entry */
						Li [nz] = k1 + kk ;
						Lx [nz] = 1 ;
						nz++ ;
						Li2 = Lx2 = GET_POINTER (LU, Lip, Lip_offset,
								Llen, Llen_offset,
								Li2_offset, Lx2_offset, kk, len) ;
						for (p = 0 ; p < len[0] ; p++)
						{
							Li [nz] = k1 + (int) Li2 [Li2_offset[0] + p] ;
							Lx [nz] = Lx2 [Lx2_offset[0] + p] ; //REAL (Lx2 [p]) ;
							nz++ ;
						}
					}
				}
			}
			Lp [n] = nz ;
			ASSERT (nz == Numeric.lnz) ;
		}

		/* ---------------------------------------------------------------------- */
		/* extract each block of U */
		/* ---------------------------------------------------------------------- */

		if (Up != null && Ui != null && Ux != null)
		{
			nz = 0 ;
			for (block = 0 ; block < nblocks ; block++)
			{
				k1 = Symbolic.R [block] ;
				k2 = Symbolic.R [block+1] ;
				nk = k2 - k1 ;
				Ukk = Numeric.Udiag ;
				int Ukk_offset = k1 ;
				if (nk == 1)
				{
					/* singleton block */
					Up [k1] = nz ;
					Ui [nz] = k1 ;
					Ux [nz] = Ukk [Ukk_offset + 0] ; //REAL (Ukk [0]) ;
					nz++ ;
				}
				else
				{
					/* non-singleton block */
					LU = Numeric.LUbx [block] ;
					Uip = Numeric.Uip ;
					int Uip_offset = k1 ;
					Ulen = Numeric.Ulen ;
					int Ulen_offset = k1 ;
					for (kk = 0 ; kk < nk ; kk++)
					{
						Up [k1+kk] = nz ;
						Ui2 = Ux2 = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui2_offset, Ux2_offset, kk, len) ;
						for (p = 0 ; p < len[0] ; p++)
						{
							Ui [nz] = k1 + (int) Ui2 [Ui2_offset[0] + p] ;
							Ux [nz] = Ux2 [Ux2_offset[0] + p] ; //REAL (Ux2 [p]) ;
							nz++ ;
						}
						/* add the diagonal entry */
						Ui [nz] = k1 + kk ;
						Ux [nz] = Ukk [Ukk_offset+ kk] ; //REAL (Ukk [kk]) ;
						nz++ ;
					}
				}
			}
			Up [n] = nz ;
			ASSERT (nz == Numeric.unz) ;
		}

		/* ---------------------------------------------------------------------- */
		/* extract the off-diagonal blocks, F */
		/* ---------------------------------------------------------------------- */

		if (Fp != null && Fi != null && Fx != null)
		{
			for (k = 0 ; k <= n ; k++)
			{
				Fp [k] = Numeric.Offp [k] ;
			}
			nz = Fp [n] ;
			for (k = 0 ; k < nz ; k++)
			{
				Fi [k] = Numeric.Offi [k] ;
			}
			for (k = 0 ; k < nz ; k++)
			{
				Fx [k] = Numeric.Offx [k] ;
				//Fx [k] = REAL (((double[]) Numeric.Offx) [k]) ;
			}
		}

		return (TRUE) ;
	}
}
