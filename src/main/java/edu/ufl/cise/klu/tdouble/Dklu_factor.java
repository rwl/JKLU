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
 * Factor the matrix, after ordering and analyzing it with KLU_analyze
 * or KLU_analyze_given.
 */
public class Dklu_factor extends Dklu_internal
{

	/**
	 *
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Ax
	 * @param Symbolic
	 * @param Numeric
	 * @param Common
	 */
	public static void factor2(int[] Ap, int[] Ai, double[] Ax,
			KLU_symbolic Symbolic, KLU_numeric Numeric, KLU_common Common)
	{
		double lsize ;
		double[] Lnz, Rs ;
		int[] P, Q, R, Pnum, Offp, Offi, Pblock, Pinv, Iwork,
			Lip, Uip, Llen, Ulen ;
		double[] Offx, X, Udiag ;
		double s ;
		double[][] LUbx ;
		int k1, k2, nk, k, block, oldcol, pend, oldrow, n, lnz, unz, p, newrow,
			nblocks, poff, nzoff, lnz_block, unz_block, scale, max_lnz_block,
			max_unz_block ;

		/* ---------------------------------------------------------------------- */
		/* initializations */
		/* ---------------------------------------------------------------------- */

		/* get the contents of the Symbolic object */
		n = Symbolic.n ;
		P = Symbolic.P ;
		Q = Symbolic.Q ;
		R = Symbolic.R ;
		Lnz = Symbolic.Lnz ;
		nblocks = Symbolic.nblocks ;
		nzoff = Symbolic.nzoff ;

		Pnum = Numeric.Pnum ;
		Offp = Numeric.Offp ;
		Offi = Numeric.Offi ;
		Offx = (double[]) Numeric.Offx ;

		Lip = Numeric.Lip ;
		Uip = Numeric.Uip ;
		Llen = Numeric.Llen ;
		Ulen = Numeric.Ulen ;
		LUbx = (double[][]) Numeric.LUbx ;
		Udiag = Numeric.Udiag ;

		Rs = Numeric.Rs ;
		Pinv = Numeric.Pinv ;
		X = (double[]) Numeric.Xwork ;              /* X is of size n */
		Iwork = Numeric.Iwork ;                    /* 5*maxblock for KLU_factor */
													/* 1*maxblock for Pblock */
		Pblock = Iwork + 5*((int) Symbolic.maxblock) ;
		Common.nrealloc = 0 ;
		scale = Common.scale ;
		max_lnz_block = 1 ;
		max_unz_block = 1 ;

		/* compute the inverse of P from symbolic analysis.  Will be updated to
		 * become the inverse of the numerical factorization when the factorization
		 * is done, for use in KLU_refactor */
		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++)
			{
				Pinv [k] = EMPTY ;
			}
		}
		for (k = 0 ; k < n ; k++)
		{
			ASSERT (P [k] >= 0 && P [k] < n) ;
			Pinv [P [k]] = k ;
		}
		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
		}

		lnz = 0 ;
		unz = 0 ;
		Common.noffdiag = 0 ;
		Offp [0] = 0 ;

		/* ---------------------------------------------------------------------- */
		/* optionally check input matrix and compute scale factors */
		/* ---------------------------------------------------------------------- */

		if (scale >= 0)
		{
			/* use Pnum as workspace. NOTE: scale factors are not yet permuted
			 * according to the final pivot row ordering, so Rs [oldrow] is the
			 * scale factor for A (oldrow,:), for the user's matrix A.  Pnum is
			 * used as workspace in KLU_scale.  When the factorization is done,
			 * the scale factors are permuted according to the final pivot row
			 * permutation, so that Rs [k] is the scale factor for the kth row of
			 * A(p,q) where p and q are the final row and column permutations. */
			Dklu_scale.klu_scale (scale, n, Ap, Ai, (double[]) Ax, Rs, Pnum, Common) ;
			if (Common.status < KLU_OK)
			{
				/* matrix is invalid */
				return ;
			}
		}

		if (!NDEBUG)
		{
			if (scale > 0)
			{
				for (k = 0 ; k < n ; k++) PRINTF ("Rs [%d] %g\n", k, Rs [k]) ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* factor each block using klu */
		/* ---------------------------------------------------------------------- */

		for (block = 0 ; block < nblocks ; block++)
		{

			/* ------------------------------------------------------------------ */
			/* the block is from rows/columns k1 to k2-1 */
			/* ------------------------------------------------------------------ */

			k1 = R [block] ;
			k2 = R [block+1] ;
			nk = k2 - k1 ;
			PRINTF ("FACTOR BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk) ;

			if (nk == 1)
			{

				/* -------------------------------------------------------------- */
				/* singleton case */
				/* -------------------------------------------------------------- */

				poff = Offp [k1] ;
				oldcol = Q [k1] ;
				pend = Ap [oldcol+1] ;
				CLEAR (s) ;

				if (scale <= 0)
				{
					/* no scaling */
					for (p = Ap [oldcol] ; p < pend ; p++)
					{
						oldrow = Ai [p] ;
						newrow = Pinv [oldrow] ;
						if (newrow < k1)
						{
							Offi [poff] = oldrow ;
							Offx [poff] = Ax [p] ;
							poff++ ;
						}
						else
						{
							ASSERT (newrow == k1) ;
							PRINTF ("singleton block %d", block) ;
							PRINT_ENTRY (Ax [p]) ;
							s = Ax [p] ;
						}
					}
				}
				else
				{
					/* row scaling.  NOTE: scale factors are not yet permuted
					 * according to the pivot row permutation, so Rs [oldrow] is
					 * used below.  When the factorization is done, the scale
					 * factors are permuted, so that Rs [newrow] will be used in
					 * klu_solve, klu_tsolve, and klu_rgrowth */
					for (p = Ap [oldcol] ; p < pend ; p++)
					{
						oldrow = Ai [p] ;
						newrow = Pinv [oldrow] ;
						if (newrow < k1)
						{
							Offi [poff] = oldrow ;
							/* Offx [poff] = Ax [p] / Rs [oldrow] ; */
							SCALE_DIV_ASSIGN (Offx [poff], Ax [p], Rs [oldrow]) ;
							poff++ ;
						}
						else
						{
							ASSERT (newrow == k1) ;
							PRINTF ("singleton block %d ", block) ;
							PRINT_ENTRY (Ax[p]) ;
							SCALE_DIV_ASSIGN (s, Ax [p], Rs [oldrow]) ;
						}
					}
				}

				Udiag [k1] = s ;

				if (IS_ZERO (s))
				{
					/* singular singleton */
					Common.status = KLU_SINGULAR ;
					Common.numerical_rank = k1 ;
					Common.singular_col = oldcol ;
					if (Common.halt_if_singular == 1)
					{
						return ;
					}
				}

				Offp [k1+1] = poff ;
				Pnum [k1] = P [k1] ;
				lnz++ ;
				unz++ ;

			}
			else
			{

				/* -------------------------------------------------------------- */
				/* construct and factorize the kth block */
				/* -------------------------------------------------------------- */

				if (Lnz [block] < 0)
				{
					/* COLAMD was used - no estimate of fill-in */
					/* use 10 times the nnz in A, plus n */
					lsize = -(Common.initmem) ;
				}
				else
				{
					lsize = Common.initmem_amd * Lnz [block] + nk ;
				}

				/* allocates 1 arrays: LUbx [block] */
				Numeric.LUsize [block] = Dklu.klu_kernel_factor (
						nk, Ap, Ai, Ax, Q,
						lsize, LUbx [block], Udiag + k1, Llen + k1, Ulen + k1,
						Lip + k1, Uip + k1, Pblock, lnz_block, unz_block,
						X, Iwork, k1, Pinv, Rs, Offp, Offi, Offx, Common) ;

				if (Common.status < KLU_OK ||
				   (Common.status == KLU_SINGULAR &&
						   Common.halt_if_singular == 1))
				{
					/* out of memory, invalid inputs, or singular */
					return ;
				}

				PRINTF ("\n----------------------- L %d:\n", block) ;
				ASSERT (Dklu_dump.klu_valid_LU (nk, TRUE, Lip+k1, Llen+k1,
						LUbx [block])) ;
				PRINTF ("\n----------------------- U %d:\n", block) ;
				ASSERT (Dklu_dump.klu_valid_LU (nk, FALSE, Uip+k1, Ulen+k1,
						LUbx [block])) ;

				/* -------------------------------------------------------------- */
				/* get statistics */
				/* -------------------------------------------------------------- */

				lnz += lnz_block ;
				unz += unz_block ;
				max_lnz_block = MAX (max_lnz_block, lnz_block) ;
				max_unz_block = MAX (max_unz_block, unz_block) ;

				if (Lnz [block] == EMPTY)
				{
					/* revise estimate for subsequent factorization */
					Lnz [block] = MAX (lnz_block, unz_block) ;
				}

				/* -------------------------------------------------------------- */
				/* combine the klu row ordering with the symbolic pre-ordering */
				/* -------------------------------------------------------------- */

				PRINTF ("Pnum, 1-based:\n") ;
				for (k = 0 ; k < nk ; k++)
				{
					ASSERT (k + k1 < n) ;
					ASSERT (Pblock [k] + k1 < n) ;
					Pnum [k + k1] = P [Pblock [k] + k1] ;
					PRINTF ("Pnum (%d + %d + 1 = %d) = %d + 1 = %d\n",
						k, k1, k+k1+1, Pnum [k+k1], Pnum [k+k1]+1) ;
				}

				/* the local pivot row permutation Pblock is no longer needed */
			}
		}
		ASSERT (nzoff == Offp [n]) ;
		PRINTF ("\n------------------- Off diagonal entries:\n") ;
		ASSERT (Dklu_dump.klu_valid (n, Offp, Offi, Offx)) ;

		Numeric.lnz = lnz ;
		Numeric.unz = unz ;
		Numeric.max_lnz_block = max_lnz_block ;
		Numeric.max_unz_block = max_unz_block ;

		/* compute the inverse of Pnum */
		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++)
			{
				Pinv [k] = EMPTY ;
			}
		}
		for (k = 0 ; k < n ; k++)
		{
			ASSERT (Pnum [k] >= 0 && Pnum [k] < n) ;
			Pinv [Pnum [k]] = k ;
		}
		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
		}

		/* permute scale factors Rs according to pivotal row order */
		if (scale > 0)
		{
			for (k = 0 ; k < n ; k++)
			{
				REAL (X [k]) = Rs [Pnum [k]] ;
			}
			for (k = 0 ; k < n ; k++)
			{
				Rs [k] = REAL (X [k]) ;
			}
		}

		PRINTF ("\n------------------- Off diagonal entries, old:\n") ;
		ASSERT (Dklu_dump.klu_valid (n, Offp, Offi, Offx)) ;

		/* apply the pivot row permutations to the off-diagonal entries */
		for (p = 0 ; p < nzoff ; p++)
		{
			ASSERT (Offi [p] >= 0 && Offi [p] < n) ;
			Offi [p] = Pinv [Offi [p]] ;
		}

		PRINTF ("\n------------------- Off diagonal entries, new:\n") ;
		ASSERT (Dklu_dump.klu_valid (n, Offp, Offi, Offx)) ;

		if (!NDEBUG)
		{
			PRINTF ("\n ############# KLU_BTF_FACTOR done, nblocks %d\n",
					nblocks);
			double ss ;
			double[] Udiag = Numeric.Udiag ;
			for (block = 0 ; block < nblocks && Common.status == KLU_OK ; block++)
			{
				k1 = R [block] ;
				k2 = R [block+1] ;
				nk = k2 - k1 ;
				PRINTF ("\n======================KLU_factor output: k1 %d k2 %d nk %d\n",k1,k2,nk) ;
				if (nk == 1)
				{
					PRINTF ("singleton  ") ;
					/* ENTRY_PRINT (singleton [block]) ; */
					ss = Udiag [k1] ;
					PRINT_ENTRY (ss) ;
				}
				else
				{
					int[] Lip, Uip, Llen, Ulen ;
					double[] LU ;
					Lip = Numeric.Lip + k1 ;
					Llen = Numeric.Llen + k1 ;
					LU = (double[]) Numeric.LUbx [block] ;
					PRINTF ("\n---- L block %d\n", block);
					ASSERT (Dklu_dump.klu_valid_LU (nk, TRUE, Lip, Llen, LU)) ;
					Uip = Numeric.Uip + k1 ;
					Ulen = Numeric.Ulen + k1 ;
					PRINTF ("\n---- U block %d\n", block) ;
					ASSERT (Dklu_dump.klu_valid_LU (nk, FALSE, Uip, Ulen, LU)) ;
				}
			}
		}
	}

	/**
	 *
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Ax
	 * @param Symbolic
	 * @param Common
	 * @return null if error, or a valid KLU_numeric object if successful
	 */
	public static KLU_numeric klu_factor(int[] Ap, int[] Ai, double[] Ax,
			KLU_symbolic Symbolic, KLU_common Common)
	{
		int n, nzoff, nblocks, maxblock, k, ok = TRUE ;
		int[] R ;
		KLU_numeric Numeric ;
		int n1, nzoff1, s, b6, n3 ;

		if (Common == null)
		{
			return (null) ;
		}
		Common.status = KLU_OK ;
		Common.numerical_rank = EMPTY ;
		Common.singular_col = EMPTY ;

		/* ---------------------------------------------------------------------- */
		/* get the contents of the Symbolic object */
		/* ---------------------------------------------------------------------- */

		/* check for a valid Symbolic object */
		if (Symbolic == null)
		{
			Common.status = KLU_INVALID ;
			return (null) ;
		}

		n = Symbolic.n ;
		nzoff = Symbolic.nzoff ;
		nblocks = Symbolic.nblocks ;
		maxblock = Symbolic.maxblock ;
		R = Symbolic.R ;
		PRINTF ("KLU_factor:  n %d nzoff %d nblocks %d maxblock %d\n",
			n, nzoff, nblocks, maxblock) ;

		/* ---------------------------------------------------------------------- */
		/* get control parameters and make sure they are in the proper range */
		/* ---------------------------------------------------------------------- */

		Common.initmem_amd = MAX (1.0, Common.initmem_amd) ;
		Common.initmem = MAX (1.0, Common.initmem) ;
		Common.tol = MIN (Common.tol, 1.0) ;
		Common.tol = MAX (0.0, Common.tol) ;
		Common.memgrow = MAX (1.0, Common.memgrow) ;

		/* ---------------------------------------------------------------------- */
		/* allocate the Numeric object  */
		/* ---------------------------------------------------------------------- */

		/* this will not cause int overflow (already checked by KLU_symbolic) */
		n1 = ((int) n) + 1 ;
		nzoff1 = ((int) nzoff) + 1 ;

		Numeric = Dklu_memory.klu_malloc (sizeof (KLU_numeric), 1, Common) ;
		if (Common.status < KLU_OK)
		{
			/* out of memory */
			Common.status = KLU_OUT_OF_MEMORY ;
			return (null) ;
		}
		Numeric.n = n ;
		Numeric.nblocks = nblocks ;
		Numeric.nzoff = nzoff ;
		Numeric.Pnum = Dklu_memory.klu_malloc (n, sizeof (Integer), Common) ;
		Numeric.Offp = Dklu_memory.klu_malloc (n1, sizeof (Integer), Common) ;
		Numeric.Offi = Dklu_memory.klu_malloc (nzoff1, sizeof (Integer), Common) ;
		Numeric.Offx = Dklu_memory.klu_malloc (nzoff1, sizeof (double), Common) ;

		Numeric.Lip  = Dklu_memory.klu_malloc (n, sizeof (Integer), Common) ;
		Numeric.Uip  = Dklu_memory.klu_malloc (n, sizeof (Integer), Common) ;
		Numeric.Llen = Dklu_memory.klu_malloc (n, sizeof (Integer), Common) ;
		Numeric.Ulen = Dklu_memory.klu_malloc (n, sizeof (Integer), Common) ;

		Numeric.LUsize = Dklu_memory.klu_malloc (nblocks, sizeof (int), Common) ;

		Numeric.LUbx = Dklu_memory.klu_malloc (nblocks, sizeof (double[]), Common) ;
		if (Numeric.LUbx != null)
		{
			for (k = 0 ; k < nblocks ; k++)
			{
				Numeric.LUbx [k] = null ;
			}
		}

		Numeric.Udiag = Dklu_memory.klu_malloc (n, sizeof (double), Common) ;

		if (Common.scale > 0)
		{
			Numeric.Rs = Dklu_memory.klu_malloc (n, sizeof (Double), Common) ;
		}
		else
		{
			/* no scaling */
			Numeric.Rs = null ;
		}

		Numeric.Pinv = Dklu_memory.klu_malloc (n, sizeof (Int), Common) ;

		/* allocate permanent workspace for factorization and solve.  Note that the
		 * solver will use an Xwork of size 4n, whereas the factorization codes use
		 * an Xwork of size n and integer space (Iwork) of size 6n. KLU_condest
		 * uses an Xwork of size 2n.  Total size is:
		 *
		 *    n*sizeof(double) + max (6*maxblock*sizeof(Int), 3*n*sizeof(double))
		 */
		s = Dklu_mult_size_t.klu_mult_size_t (n, sizeof (double), ok) ;
		n3 = Dklu_mult_size_t.klu_mult_size_t (n, 3 * sizeof (double), ok) ;
		b6 = Dklu_mult_size_t.klu_mult_size_t (maxblock, 6 * sizeof (Int), ok) ;
		Numeric.worksize = Dklu_add_size_t.klu_add_size_t (s, MAX (n3, b6), ok) ;
		Numeric.Work = Dklu_memory.klu_malloc (Numeric.worksize, 1, Common) ;
		Numeric.Xwork = Numeric.Work ;
		Numeric.Iwork = (Int[]) ((double[]) Numeric.Xwork + n) ;
		if (!ok || Common.status < KLU_OK)
		{
			/* out of memory or problem too large */
			Common.status = ok == 1 ? KLU_OUT_OF_MEMORY :
					KLU_TOO_LARGE ;
			Dklu_free_numeric.klu_free_numeric (Numeric, Common) ;
			return (null) ;
		}

		/* ---------------------------------------------------------------------- */
		/* factorize the blocks */
		/* ---------------------------------------------------------------------- */

		factor2 (Ap, Ai, (double[]) Ax, Symbolic, Numeric, Common) ;

		/* ---------------------------------------------------------------------- */
		/* return or free the Numeric object */
		/* ---------------------------------------------------------------------- */

		if (Common.status < KLU_OK)
		{
			/* out of memory or inputs invalid */
			Dklu_free_numeric.klu_free_numeric (Numeric, Common) ;
		}
		else if (Common.status == KLU_SINGULAR)
		{
			if (Common.halt_if_singular == 1)
			{
				/* Matrix is singular, and the Numeric object is only partially
				 * defined because we halted early.  This is the default case for
				 * a singular matrix. */
				Dklu_free_numeric.klu_free_numeric (Numeric, Common) ;
			}
		}
		else if (Common.status == KLU_OK)
		{
			/* successful non-singular factorization */
			Common.numerical_rank = n ;
			Common.singular_col = n ;
		}
		return (Numeric) ;
	}

}
