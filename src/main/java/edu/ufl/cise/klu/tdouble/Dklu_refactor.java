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
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;

import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_malloc_dbl;
import static edu.ufl.cise.klu.tdouble.Dklu_scale.klu_scale;
import static edu.ufl.cise.klu.tdouble.Dklu_dump.klu_valid;
import static edu.ufl.cise.klu.tdouble.Dklu_dump.klu_valid_LU;

/**
 * Factor the matrix, after ordering and analyzing it with KLU_analyze, and
 * factoring it once with KLU_factor.  This routine cannot do any numerical
 * pivoting.  The pattern of the input matrix (Ap, Ai) must be identical to
 * the pattern given to KLU_factor.
 */
public class Dklu_refactor extends Dklu_internal {

	/**
	 *
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Ax
	 * @param Symbolic
	 * @param Numeric
	 * @param Common
	 * @return true if successful, false otherwise
	 */
	public static int klu_refactor(int[] Ap, int[] Ai, double[] Ax,
			KLU_symbolic Symbolic, KLU_numeric Numeric, KLU_common  Common)
	{
		double ukk, ujk, s ;
		double[] Offx, Lx, Ux, X, Az, Udiag ;
		double[] Rs ;
		int[] Q, R, Pnum, Offp, Offi, Pinv, Lip, Uip, Llen, Ulen ;
		/*int[]*/double[] Ui, Li ;
		double[][] LUbx ;
		double[] LU ;
		int k1, k2, nk, k, block, oldcol, pend, oldrow, n, p, newrow, scale,
			nblocks, poff, i, j, up, maxblock, nzoff ;
		int[] ulen = new int[1] ;
		int[] Ui_offset = new int[1] ;
		int[] Ux_offset = new int[1] ;
		int[] llen = new int[1] ;
		int[] Li_offset = new int[1] ;
		int[] Lx_offset = new int[1] ;

		/* ---------------------------------------------------------------------- */
		/* check inputs */
		/* ---------------------------------------------------------------------- */

		if (Common == null)
		{
			return (FALSE) ;
		}
		Common.status = KLU_OK ;

		if (Numeric == null)
		{
			/* invalid Numeric object */
			Common.status = KLU_INVALID ;
			return (FALSE) ;
		}

		Common.numerical_rank = EMPTY ;
		Common.singular_col = EMPTY ;

		Az = (double[]) Ax ;

		/* ---------------------------------------------------------------------- */
		/* get the contents of the Symbolic object */
		/* ---------------------------------------------------------------------- */

		n = Symbolic.n ;
		Q = Symbolic.Q ;
		R = Symbolic.R ;
		nblocks = Symbolic.nblocks ;
		maxblock = Symbolic.maxblock ;

		/* ---------------------------------------------------------------------- */
		/* get the contents of the Numeric object */
		/* ---------------------------------------------------------------------- */

		Pnum = Numeric.Pnum ;
		Offp = Numeric.Offp ;
		Offi = Numeric.Offi ;
		Offx = Numeric.Offx ;

		LUbx = Numeric.LUbx ;

		scale = Common.scale ;
		if (scale > 0)
		{
			/* factorization was not scaled, but refactorization is scaled */
			if (Numeric.Rs == null)
			{
				Numeric.Rs = klu_malloc_dbl (n, Common) ;
				if (Common.status < KLU_OK)
				{
					Common.status = KLU_OUT_OF_MEMORY ;
					return (FALSE) ;
				}
			}
		}
		else
		{
			/* no scaling for refactorization; ensure Numeric.Rs is freed.  This
			 * does nothing if Numeric.Rs is already null. */
			//Numeric.Rs = KLU_free (Numeric.Rs, n, sizeof (double), Common) ;
			Numeric.Rs = null ;
		}
		Rs = Numeric.Rs ;

		Pinv = Numeric.Pinv ;
		X = Numeric.Xwork ;
		Common.nrealloc = 0 ;
		Udiag = Numeric.Udiag ;
		nzoff = Symbolic.nzoff ;

		/* ---------------------------------------------------------------------- */
		/* check the input matrix compute the row scale factors, Rs */
		/* ---------------------------------------------------------------------- */

		/* do no scale, or check the input matrix, if scale < 0 */
		if (scale >= 0)
		{
			/* check for out-of-range indices, but do not check for duplicates */
			if (klu_scale (scale, n, Ap, Ai, Ax, Rs, null, Common) == 0)
			{
				return (FALSE) ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* clear workspace X */
		/* ---------------------------------------------------------------------- */

		for (k = 0 ; k < maxblock ; k++)
		{
			/* X [k] = 0 ; */
			CLEAR (X, k) ;
		}

		poff = 0 ;

		/* ---------------------------------------------------------------------- */
		/* factor each block */
		/* ---------------------------------------------------------------------- */

		if (scale <= 0)
		{

			/* ------------------------------------------------------------------ */
			/* no scaling */
			/* ------------------------------------------------------------------ */

			for (block = 0 ; block < nblocks ; block++)
			{

				/* -------------------------------------------------------------- */
				/* the block is from rows/columns k1 to k2-1 */
				/* -------------------------------------------------------------- */

				k1 = R [block] ;
				k2 = R [block+1] ;
				nk = k2 - k1 ;

				if (nk == 1)
				{

					/* ---------------------------------------------------------- */
					/* singleton case */
					/* ---------------------------------------------------------- */

					oldcol = Q [k1] ;
					pend = Ap [oldcol+1] ;
					s = 0 ; //CLEAR (s) ;
					for (p = Ap [oldcol] ; p < pend ; p++)
					{
						newrow = Pinv [Ai [p]] - k1 ;
						if (newrow < 0 && poff < nzoff)
						{
							/* entry in off-diagonal block */
							Offx [poff] = Az [p] ;
							poff++ ;
						}
						else
						{
							/* singleton */
							s = Az [p] ;
						}
					}
					Udiag [k1] = s ;

				}
				else
				{

					/* ---------------------------------------------------------- */
					/* construct and factor the kth block */
					/* ---------------------------------------------------------- */

					Lip  = Numeric.Lip ;
					int Lip_offset = k1 ;
					Llen = Numeric.Llen ;
					int Llen_offset = k1 ;
					Uip  = Numeric.Uip ;
					int Uip_offset = k1 ;
					Ulen = Numeric.Ulen ;
					int Ulen_offset = k1 ;
					LU = LUbx [block] ;

					for (k = 0 ; k < nk ; k++)
					{

						/* ------------------------------------------------------ */
						/* scatter kth column of the block into workspace X */
						/* ------------------------------------------------------ */

						oldcol = Q [k+k1] ;
						pend = Ap [oldcol+1] ;
						for (p = Ap [oldcol] ; p < pend ; p++)
						{
							newrow = Pinv [Ai [p]] - k1 ;
							if (newrow < 0 && poff < nzoff)
							{
								/* entry in off-diagonal block */
								Offx [poff] = Az [p] ;
								poff++ ;
							}
							else
							{
								/* (newrow,k) is an entry in the block */
								X [newrow] = Az [p] ;
							}
						}

						/* ------------------------------------------------------ */
						/* compute kth column of U, and update kth column of A */
						/* ------------------------------------------------------ */

						Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
								Ui_offset, Ux_offset, k, ulen) ;
						for (up = 0 ; up < ulen[0] ; up++)
						{
							j = (int) Ui [Ui_offset[0] + up] ;
							ujk = X [j] ;
							/* X [j] = 0 ; */
							CLEAR (X, j) ;
							Ux [Ux_offset[0] + up] = ujk ;
							Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset,
									Li_offset, Lx_offset, j, llen) ;
							for (p = 0 ; p < llen[0] ; p++)
							{
								//MULT_SUB (X [Li [p]], Lx [p], ujk) ;
								X [(int) Li [Li_offset[0] + p]] -= Lx [Lx_offset[0] + p] * ujk ;
							}
						}
						/* get the diagonal entry of U */
						ukk = X [k] ;
						/* X [k] = 0 ; */
						CLEAR (X, k) ;
						/* singular case */
						if (IS_ZERO (ukk))
						{
							/* matrix is numerically singular */
							Common.status = KLU_SINGULAR ;
							if (Common.numerical_rank == EMPTY)
							{
								Common.numerical_rank = k+k1 ;
								Common.singular_col = Q [k+k1] ;
							}
							if (Common.halt_if_singular != 0)
							{
								/* do not continue the factorization */
								return (FALSE) ;
							}
						}
						Udiag [k+k1] = ukk ;
						/* gather and divide by pivot to get kth column of L */
						Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset,
								Li_offset, Lx_offset, k, llen) ;
						for (p = 0 ; p < llen[0] ; p++)
						{
							i = (int) Li [Li_offset[0] + p] ;
							//DIV (Lx [p], X [i], ukk) ;
							Lx [Lx_offset[0] + p] = X [i] / ukk ;
							CLEAR (X, i) ;
						}

					}
				}
			}

		}
		else
		{

			/* ------------------------------------------------------------------ */
			/* scaling */
			/* ------------------------------------------------------------------ */

			for (block = 0 ; block < nblocks ; block++)
			{

				/* -------------------------------------------------------------- */
				/* the block is from rows/columns k1 to k2-1 */
				/* -------------------------------------------------------------- */

				k1 = R [block] ;
				k2 = R [block+1] ;
				nk = k2 - k1 ;

				if (nk == 1)
				{

					/* ---------------------------------------------------------- */
					/* singleton case */
					/* ---------------------------------------------------------- */

					oldcol = Q [k1] ;
					pend = Ap [oldcol+1] ;
					s = 0 ; //CLEAR (s) ;
					for (p = Ap [oldcol] ; p < pend ; p++)
					{
						oldrow = Ai [p] ;
						newrow = Pinv [oldrow] - k1 ;
						if (newrow < 0 && poff < nzoff)
						{
							/* entry in off-diagonal block */
							Offx [poff] = Az [p] / Rs [oldrow] ;
							//SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]) ;
							poff++ ;
						}
						else
						{
							/* singleton */
							s = Az [p] / Rs [oldrow] ;
							//SCALE_DIV_ASSIGN (s, Az [p], Rs [oldrow]) ;
						}
					}
					Udiag [k1] = s ;

				}
				else
				{

					/* ---------------------------------------------------------- */
					/* construct and factor the kth block */
					/* ---------------------------------------------------------- */

					Lip  = Numeric.Lip ;
					int Lip_offset = k1 ;
					Llen = Numeric.Llen ;
					int Llen_offset = k1 ;
					Uip  = Numeric.Uip ;
					int Uip_offset = k1 ;
					Ulen = Numeric.Ulen ;
					int Ulen_offset = k1 ;
					LU = LUbx [block] ;

					for (k = 0 ; k < nk ; k++)
					{

						/* ------------------------------------------------------ */
						/* scatter kth column of the block into workspace X */
						/* ------------------------------------------------------ */

						oldcol = Q [k+k1] ;
						pend = Ap [oldcol+1] ;
						for (p = Ap [oldcol] ; p < pend ; p++)
						{
							oldrow = Ai [p] ;
							newrow = Pinv [oldrow] - k1 ;
							if (newrow < 0 && poff < nzoff)
							{
								/* entry in off-diagonal part */
								//SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]);
								Offx [poff] = Az [p] / Rs [oldrow] ;
								poff++ ;
							}
							else
							{
								/* (newrow,k) is an entry in the block */
								//SCALE_DIV_ASSIGN (X [newrow], Az [p], Rs [oldrow]) ;
								X [newrow] = Az [p] / Rs [oldrow] ;
							}
						}

						/* ------------------------------------------------------ */
						/* compute kth column of U, and update kth column of A */
						/* ------------------------------------------------------ */

						Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
								Ui_offset, Ux_offset, k, ulen) ;
						for (up = 0 ; up < ulen[0] ; up++)
						{
							j = (int) Ui [Ui_offset[0] + up] ;
							ujk = X [j] ;
							/* X [j] = 0 ; */
							CLEAR (X, j) ;
							Ux [Ux_offset[0] + up] = ujk ;
							Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset,
									Li_offset, Lx_offset, j, llen) ;
							for (p = 0 ; p < llen[0] ; p++)
							{
								//MULT_SUB (X [Li [p]], Lx [p], ujk) ;
								X [(int) Li [Li_offset[0] + p]] -= Lx [Lx_offset[0] + p] * ujk ;
							}
						}
						/* get the diagonal entry of U */
						ukk = X [k] ;
						/* X [k] = 0 ; */
						CLEAR (X, k) ;
						/* singular case */
						if (IS_ZERO (ukk))
						{
							/* matrix is numerically singular */
							Common.status = KLU_SINGULAR ;
							if (Common.numerical_rank == EMPTY)
							{
								Common.numerical_rank = k+k1 ;
								Common.singular_col = Q [k+k1] ;
							}
							if (Common.halt_if_singular != 0)
							{
								/* do not continue the factorization */
								return (FALSE) ;
							}
						}
						Udiag [k+k1] = ukk ;
						/* gather and divide by pivot to get kth column of L */
						Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset,
								Li_offset, Lx_offset, k, llen) ;
						for (p = 0 ; p < llen[0] ; p++)
						{
							i = (int) Li [Li_offset[0] + p] ;
							//DIV (Lx [p], X [i], ukk) ;
							Lx [Lx_offset[0] + p] = X [i] / ukk ;
							CLEAR (X, i) ;
						}
					}
				}
			}
		}

		/* ---------------------------------------------------------------------- */
		/* permute scale factors Rs according to pivotal row order */
		/* ---------------------------------------------------------------------- */

		if (scale > 0)
		{
			for (k = 0 ; k < n ; k++)
			{
				X [k] = Rs [Pnum [k]] ;
				//REAL (X [k]) = Rs [Pnum [k]] ;
			}
			for (k = 0 ; k < n ; k++)
			{
				Rs [k] = X [k] ;
				//Rs [k] = REAL (X [k]) ;
			}
		}

		if (!NDEBUG)
		{
			ASSERT (Offp [n] == poff) ;
			ASSERT (Symbolic.nzoff == poff) ;
			PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
			if (!NDEBUG) ASSERT (klu_valid (n, Offp, Offi, Offx)) ;
			if (Common.status == KLU_OK)
			{
				PRINTF ("\n ########### KLU_BTF_REFACTOR done, nblocks %d\n",
						nblocks);
				for (block = 0 ; block < nblocks ; block++)
				{
					k1 = R [block] ;
					k2 = R [block+1] ;
					nk = k2 - k1 ;
					PRINTF (
						"\n================KLU_refactor output: k1 %d k2 %d nk %d\n",
						k1, k2, nk) ;
					if (nk == 1)
					{
						PRINTF ("singleton  ") ;
						PRINT_ENTRY (Udiag [k1]) ;
					}
					else
					{
						Lip = Numeric.Lip ;
						int Lip_offset = k1 ;
						Llen = Numeric.Llen ;
						int Llen_offset = k1 ;
						LU = (double[]) Numeric.LUbx [block] ;
						PRINTF ("\n---- L block %d\n", block) ;
						if (!NDEBUG) ASSERT (klu_valid_LU (nk, TRUE, Lip, Lip_offset, Llen, Llen_offset, LU)) ;
						Uip = Numeric.Uip ;
						int Uip_offset = k1 ;
						Ulen = Numeric.Ulen ;
						int Ulen_offset = k1 ;
						PRINTF ("\n---- U block %d\n", block) ;
						if (!NDEBUG) ASSERT (klu_valid_LU (nk, FALSE, Uip, Uip_offset, Ulen, Ulen_offset, LU)) ;
					}
				}
			}
		}

		return (TRUE) ;
	}

}
