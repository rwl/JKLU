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
import edu.ufl.cise.klu.common.KLU_symbolic;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze_given.klu_analyze_given;
import static edu.ufl.cise.klu.tdouble.Dklu_analyze_given.klu_alloc_symbolic;
import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_malloc_int;
import static edu.ufl.cise.klu.tdouble.Dklu_dump.klu_valid;

import static edu.ufl.cise.amd.tdouble.Damd_order.amd_order;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_INFO;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_OK;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_OUT_OF_MEMORY;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_MEMORY;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_LNZ;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_NMULTSUBS_LU;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_NDIV;
import static edu.ufl.cise.amd.tdouble.Damd.AMD_SYMMETRY;

import static edu.ufl.cise.colamd.tdouble.Dcolamd.colamd;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_recommended;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_STATS;

import static edu.ufl.cise.btf.tdouble.Dbtf_order.btf_order;
import static edu.ufl.cise.btf.tdouble.Dbtf.BTF_UNFLIP;

/**
 * Orders and analyzes a matrix.
 *
 * Order the matrix using BTF (or not), and then AMD, COLAMD, the natural
 * ordering, or the user-provided-function on the blocks.  Does not support
 * using a given ordering (use klu_analyze_given for that case).
 */
public class Dklu_analyze extends Dklu_internal
{

	/**
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param nblocks # of blocks
	 * @param Pbtf BTF row permutation
	 * @param Qbtf BTF col permutation
	 * @param R size n+1, but only Rbtf [0..nblocks] is used
	 * @param ordering what ordering to use (0, 1, or 3 for this routine)
	 * @param P size n
	 * @param Q size n
	 * @param Lnz size n, but only Lnz [0..nblocks-1] is used
	 * @param Pblk size maxblock
	 * @param Cp size maxblock+1
	 * @param Ci size MAX (nz+1, Cilen)
	 * @param Cilen nz+1, or COLAMD_recommend(nz,n,n) for COLAMD
	 * @param Pinv size maxblock
	 * @param Symbolic
	 * @param Common
	 * @return KLU_OK or < 0 if error
	 */
	public static int analyze_worker(int n, int[] Ap, int[] Ai, int nblocks,
			int[] Pbtf, int[] Qbtf, int[] R, int ordering, int[] P, int[] Q,
			double[] Lnz, int[] Pblk, int[] Cp, int[] Ci, int Cilen,
			int[] Pinv, KLU_symbolic Symbolic, KLU_common Common)
	{
		double[] amd_Info = new double[AMD_INFO] ;
		double lnz, lnz1, flops, flops1 ;
		int k1, k2, nk, k, block, oldcol, pend, newcol, result, pc, p, newrow,
			maxnz, nzoff, ok, err = KLU_INVALID ;
		int[] cstats = new int[COLAMD_STATS];

		/* ---------------------------------------------------------------------- */
		/* initializations */
		/* ---------------------------------------------------------------------- */

		/* compute the inverse of Pbtf */
		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++)
			{
				P [k] = EMPTY ;
				Q [k] = EMPTY ;
				Pinv [k] = EMPTY ;
			}
		}
		for (k = 0 ; k < n ; k++)
		{
			ASSERT (Pbtf [k] >= 0 && Pbtf [k] < n) ;
			Pinv [Pbtf [k]] = k ;
		}
		if (!NDEBUG) {
			for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
		}
		nzoff = 0 ;
		lnz = 0 ;
		maxnz = 0 ;
		flops = 0 ;
		Symbolic.symmetry = EMPTY ;        /* only computed by AMD */

		/* ---------------------------------------------------------------------- */
		/* order each block */
		/* ---------------------------------------------------------------------- */

		for (block = 0 ; block < nblocks ; block++)
		{

			/* ------------------------------------------------------------------ */
			/* the block is from rows/columns k1 to k2-1 */
			/* ------------------------------------------------------------------ */

			k1 = R [block] ;
			k2 = R [block+1] ;
			nk = k2 - k1 ;
			PRINTF ("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk) ;

			/* ------------------------------------------------------------------ */
			/* construct the kth block, C */
			/* ------------------------------------------------------------------ */

			Lnz [block] = EMPTY ;
			pc = 0 ;
			for (k = k1 ; k < k2 ; k++)
			{
				newcol = k-k1 ;
				Cp [newcol] = pc ;
				oldcol = Qbtf [k] ;
				pend = Ap [oldcol+1] ;
				for (p = Ap [oldcol] ; p < pend ; p++)
				{
					newrow = Pinv [Ai [p]] ;
					if (newrow < k1)
					{
						nzoff++ ;
					}
					else
					{
						/* (newrow, newcol) is an entry in the block */
						ASSERT (newrow < k2) ;
						newrow -= k1 ;
						Ci [pc++] = newrow ;
					}
				}
			}
			Cp [nk] = pc ;
			maxnz = MAX (maxnz, pc) ;
			if (!NDEBUG) ASSERT (klu_valid (nk, Cp, Ci, null)) ;

			/* ------------------------------------------------------------------ */
			/* order the block C */
			/* ------------------------------------------------------------------ */

			if (nk <= 3)
			{

				/* -------------------------------------------------------------- */
				/* use natural ordering for tiny blocks (3-by-3 or less) */
				/* -------------------------------------------------------------- */

				for (k = 0 ; k < nk ; k++)
				{
					Pblk [k] = k ;
				}
				lnz1 = nk * (nk + 1) / 2 ;
				flops1 = nk * (nk - 1) / 2 + (nk-1)*nk*(2*nk-1) / 6 ;
				ok = TRUE ;

			}
			else if (ordering == 0)
			{

				/* -------------------------------------------------------------- */
				/* order the block with AMD (C+C') */
				/* -------------------------------------------------------------- */

				result = amd_order (nk, Cp, Ci, Pblk, null, amd_Info) ;
				ok = (result >= AMD_OK) ? 1 : 0;
				if (result == AMD_OUT_OF_MEMORY)
				{
					err = KLU_OUT_OF_MEMORY ;
				}

				/* account for memory usage in AMD */
				Common.mempeak = MAX (Common.mempeak,
					Common.memusage + (long) amd_Info [AMD_MEMORY]) ;

				/* get the ordering statistics from AMD */
				lnz1 = (int) (amd_Info [AMD_LNZ]) + nk ;
				flops1 = 2 * amd_Info [AMD_NMULTSUBS_LU] + amd_Info [AMD_NDIV] ;
				if (pc == maxnz)
				{
					/* get the symmetry of the biggest block */
					Symbolic.symmetry = amd_Info [AMD_SYMMETRY] ;
				}

			}
			else if (ordering == 1)
			{

				/* -------------------------------------------------------------- */
				/* order the block with COLAMD (C) */
				/* -------------------------------------------------------------- */

				/* order (and destroy) Ci, returning column permutation in Cp.
				 * COLAMD "cannot" fail since the matrix has already been checked,
				 * and Ci allocated. */

				ok = colamd (nk, nk, Cilen, Ci, Cp, null, cstats) ;
				lnz1 = EMPTY ;
				flops1 = EMPTY ;

				/* copy the permutation from Cp to Pblk */
				for (k = 0 ; k < nk ; k++)
				{
					Pblk [k] = Cp [k] ;
				}

			}
			else
			{

				/* -------------------------------------------------------------- */
				/* pass the block to the user-provided ordering function */
				/* -------------------------------------------------------------- */

				lnz1 = Common.user_order.order(nk, Cp, Ci, Pblk, Common) ;
				flops1 = EMPTY ;
				ok = (lnz1 != 0) ? 1 : 0 ;
			}

			if (ok != 1)
			{
				return (err) ;  /* ordering method failed */
			}

			/* ------------------------------------------------------------------ */
			/* keep track of nnz(L) and flops statistics */
			/* ------------------------------------------------------------------ */

			Lnz [block] = lnz1 ;
			lnz = (lnz == EMPTY || lnz1 == EMPTY) ? EMPTY : (lnz + lnz1) ;
			flops = (flops == EMPTY || flops1 == EMPTY) ? EMPTY : (flops + flops1) ;

			/* ------------------------------------------------------------------ */
			/* combine the preordering with the BTF ordering */
			/* ------------------------------------------------------------------ */

			PRINTF ("Pblk, 1-based:\n") ;
			for (k = 0 ; k < nk ; k++)
			{
				ASSERT (k + k1 < n) ;
				ASSERT (Pblk [k] + k1 < n) ;
				Q [k + k1] = Qbtf [Pblk [k] + k1] ;
			}
			for (k = 0 ; k < nk ; k++)
			{
				ASSERT (k + k1 < n) ;
				ASSERT (Pblk [k] + k1 < n) ;
				P [k + k1] = Pbtf [Pblk [k] + k1] ;
			}
		}

		PRINTF ("nzoff %d  Ap[n] %d\n", nzoff, Ap [n]) ;
		ASSERT (nzoff >= 0 && nzoff <= Ap [n]) ;

		/* return estimates of # of nonzeros in L including diagonal */
		Symbolic.lnz = lnz ;           /* EMPTY if COLAMD used */
		Symbolic.unz = lnz ;
		Symbolic.nzoff = nzoff ;
		Symbolic.est_flops = flops ;   /* EMPTY if COLAMD or user-ordering used */
		return (KLU_OK) ;
	}

	/**
	 * Orders the matrix with or with BTF, then orders each block with AMD, COLAMD,
	 * or the user ordering function.  Does not handle the natural or given
	 * ordering cases.
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Common
	 * @return null if error, or a valid KLU_symbolic object if successful
	 */
	public static KLU_symbolic order_and_analyze(int n, int[] Ap, int[] Ai, KLU_common Common)
	{
		double[] work = new double [1] ;
		KLU_symbolic Symbolic ;
		double[] Lnz ;
		int[] structural_rank = new int [1] ;
		int[] Qbtf, Cp, Ci, Pinv, Pblk, Pbtf, P, Q, R ;
		int nblocks, nz, block, maxblock, k1, k2, nk, do_btf, ordering, k, Cilen ;

		/* ---------------------------------------------------------------------- */
		/* allocate the Symbolic object, and check input matrix */
		/* ---------------------------------------------------------------------- */

		Symbolic = klu_alloc_symbolic (n, Ap, Ai, Common) ;
		if (Symbolic == null)
		{
			return (null) ;
		}
		P = Symbolic.P ;
		Q = Symbolic.Q ;
		R = Symbolic.R ;
		Lnz = Symbolic.Lnz ;
		nz = Symbolic.nz ;

		ordering = Common.ordering ;
		if (ordering == 1)
		{
			/* COLAMD */
			Cilen = COLAMD_recommended (nz, n, n) ;
		}
		else if (ordering == 0 || (ordering == 3 && Common.user_order != null))
		{
			/* AMD or user ordering function */
			Cilen = nz+1 ;
		}
		else
		{
			/* invalid ordering */
			Common.status = KLU_INVALID ;
			//klu_free_symbolic (Symbolic, Common) ;
			Symbolic = null;
			return (null) ;
		}

		/* AMD memory management routines */
		//amd_malloc  = Common.malloc_memory ;
		//amd_free    = Common.free_memory ;
		//amd_calloc  = Common.calloc_memory ;
		//amd_realloc = Common.realloc_memory ;

		/* ---------------------------------------------------------------------- */
		/* allocate workspace for BTF permutation */
		/* ---------------------------------------------------------------------- */

		Pbtf = klu_malloc_int (n, Common) ;
		Qbtf = klu_malloc_int (n, Common) ;
		if (Common.status < KLU_OK)
		{
			//KLU_free (Pbtf, n, sizeof (int), Common) ;
			Pbtf = null;
			//KLU_free (Qbtf, n, sizeof (int), Common) ;
			Qbtf = null;
			//klu_free_symbolic (Symbolic, Common) ;
			Symbolic = null;
			return (null) ;
		}

		/* ---------------------------------------------------------------------- */
		/* get the common parameters for BTF and ordering method */
		/* ---------------------------------------------------------------------- */

		do_btf = Common.btf ;
		do_btf = (do_btf != 0) ? TRUE : FALSE ;
		Symbolic.ordering = ordering ;
		Symbolic.do_btf = do_btf ;
		Symbolic.structural_rank = EMPTY ;

		/* ---------------------------------------------------------------------- */
		/* find the block triangular form (if requested) */
		/* ---------------------------------------------------------------------- */

		Common.work = 0 ;

		if (do_btf != 0)
		{
			//Work = klu_malloc_int (5*n, Common) ;
			if (Common.status < KLU_OK)
			{
				/* out of memory */
				//klu_free (Pbtf, n, sizeof (int), Common) ;
				Pbtf = null;
				//klu_free (Qbtf, n, sizeof (int), Common) ;
				Qbtf = null;
				//klu_free_symbolic (Symbolic, Common) ;
				Symbolic = null;
				return (null) ;
			}

			nblocks = btf_order (n, Ap, Ai, Common.maxwork, work, Pbtf, Qbtf, R,
					structural_rank) ;
			Symbolic.structural_rank = structural_rank[0] ;
			Common.structural_rank = Symbolic.structural_rank ;
			Common.work += work[0] ;

			//klu_free (Work, 5*n, sizeof (int), Common) ;

			/* unflip Qbtf if the matrix does not have full structural rank */
			if (Symbolic.structural_rank < n)
			{
				for (k = 0 ; k < n ; k++)
				{
					Qbtf [k] = BTF_UNFLIP (Qbtf [k]) ;
				}
			}

			/* find the size of the largest block */
			maxblock = 1 ;
			for (block = 0 ; block < nblocks ; block++)
			{
				k1 = R [block] ;
				k2 = R [block+1] ;
				nk = k2 - k1 ;
				PRINTF ("block %d size %d\n", block, nk) ;
				maxblock = MAX (maxblock, nk) ;
			}
		}
		else
		{
			/* BTF not requested */
			nblocks = 1 ;
			maxblock = n ;
			R [0] = 0 ;
			R [1] = n ;
			for (k = 0 ; k < n ; k++)
			{
				Pbtf [k] = k ;
				Qbtf [k] = k ;
			}
		}

		Symbolic.nblocks = nblocks ;

		PRINTF ("maxblock size %d\n", maxblock) ;
		Symbolic.maxblock = maxblock ;

		/* ---------------------------------------------------------------------- */
		/* allocate more workspace, for analyze_worker */
		/* ---------------------------------------------------------------------- */

		Pblk = klu_malloc_int (maxblock, Common) ;
		Cp   = klu_malloc_int (maxblock + 1, Common) ;
		Ci   = klu_malloc_int (MAX (Cilen, nz+1), Common) ;
		Pinv = klu_malloc_int (n, Common) ;

		/* ---------------------------------------------------------------------- */
		/* order each block of the BTF ordering, and a fill-reducing ordering */
		/* ---------------------------------------------------------------------- */

		if (Common.status == KLU_OK)
		{
			PRINTF (("calling analyze_worker\n")) ;
			Common.status = analyze_worker (n, Ap, Ai, nblocks, Pbtf, Qbtf, R,
				ordering, P, Q, Lnz, Pblk, Cp, Ci, Cilen, Pinv, Symbolic, Common) ;
			PRINTF ("analyze_worker done\n") ;
		}

		/* ---------------------------------------------------------------------- */
		/* free all workspace */
		/* ---------------------------------------------------------------------- */

		//klu_free (Pblk, maxblock, sizeof (int), Common) ;
		//klu_free (Cp, maxblock+1, sizeof (int), Common) ;
		//klu_free (Ci, MAX (Cilen, nz+1), sizeof (int), Common) ;
		//klu_free (Pinv, n, sizeof (int), Common) ;
		//klu_free (Pbtf, n, sizeof (int), Common) ;
		//klu_free (Qbtf, n, sizeof (int), Common) ;
		Pblk = Cp = Ci = Pinv = Pbtf = Qbtf = null;

		/* ---------------------------------------------------------------------- */
		/* return the symbolic object */
		/* ---------------------------------------------------------------------- */

		if (Common.status < KLU_OK)
		{
			//klu_free_symbolic (Symbolic, Common) ;
			Symbolic = null;
		}
		return (Symbolic) ;
	}

	/**
	 * Order the matrix with BTF (or not), then order each block with AMD,
	 * COLAMD, a natural ordering, or with a user-provided ordering function.
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Common
	 * @return null if error, or a valid KLU_symbolic object if successful
	 */
	public static KLU_symbolic klu_analyze(int n, int[] Ap, int[] Ai, KLU_common Common)
	{
		/* ---------------------------------------------------------------------- */
		/* get the control parameters for BTF and ordering method */
		/* ---------------------------------------------------------------------- */

		if (Common == null)
		{
			return (null) ;
		}
		Common.status = KLU_OK ;
		Common.structural_rank = EMPTY ;

		/* ---------------------------------------------------------------------- */
		/* order and analyze */
		/* ---------------------------------------------------------------------- */

		if (Common.ordering == 2)
		{
			/* natural ordering */
			return (klu_analyze_given (n, Ap, Ai, null, null, Common)) ;
		}
		else
		{
			/* order with P and Q */
			return (order_and_analyze (n, Ap, Ai, Common)) ;
		}
	}

}
