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

import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_realloc;

/**
 * Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 */
public class Dklu_kernel extends Dklu_internal {

	/**
	 * Does a depth-first-search, starting at node j.
	 *
	 * @param j node at which to start the DFS
	 * @param k mark value, for the Flag array
	 * @param Pinv Pinv[i] = k if row i is kth pivot row, or EMPTY if
	 * row i is not yet pivotal.
	 * @param Llen size n, Llen[k] = # nonzeros in column k of L
	 * @param Lip size n, Lip[k] is position in LU of column k of L
	 * @param Stack size n
	 * @param Flag Flag[i] == k means i is marked
	 * @param Lpend for symmetric pruning
	 * @param top top of stack on input
	 * @param LU
	 * @param Lik Li row index array of the kth column
	 * @param plength
	 * @param Ap_pos keeps track of position in adj list during DFS
	 * @return
	 */
	public static int dfs(int j, int k, int[] Pinv, int[] Llen, int[] Lip,
			int[] Stack, int[] Flag, int[] Lpend, int top, double LU,
			int[] Lik, int plength, int[] Ap_pos)
	{
		int i, pos, jnew, head, l_length;
		int[] Li;

		l_length = plength ;

		head = 0 ;
		Stack [0] = j ;
		ASSERT (Flag [j] != k) ;

		while (head >= 0)
		{
			j = Stack [head] ;
			jnew = Pinv [j] ;
			ASSERT (jnew >= 0 && jnew < k) ;        /* j is pivotal */

			if (Flag [j] != k)          /* a node is not yet visited */
			{
				/* first time that j has been visited */
				Flag [j] = k ;
				PRINTF ("[ start dfs at %d : new %d\n", j, jnew) ;
				/* set Ap_pos [head] to one past the last entry in col j to scan */
				Ap_pos [head] =
					(Lpend [jnew] == EMPTY) ?  Llen [jnew] : Lpend [jnew] ;
			}

			/* add the adjacent nodes to the recursive stack by iterating through
			 * until finding another non-visited pivotal node */
			Li = (int[]) (LU + Lip [jnew]) ;
			for (pos = --Ap_pos [head] ; pos >= 0 ; --pos)
			{
				i = Li [pos] ;
				if (Flag [i] != k)
				{
					/* node i is not yet visited */
					if (Pinv [i] >= 0)
					{
						/* keep track of where we left off in the scan of the
						 * adjacency list of node j so we can restart j where we
						 * left off. */
						Ap_pos [head] = pos ;

						/* node i is pivotal; push it onto the recursive stack
						 * and immediately break so we can recurse on node i. */
						Stack [++head] = i ;
						break ;
					}
					else
					{
						/* node i is not pivotal (no outgoing edges). */
						/* Flag as visited and store directly into L,
						 * and continue with current node j. */
						Flag [i] = k ;
						Lik [l_length] = i ;
						l_length++ ;
					}
				}
			}

			if (pos == -1)
			{
				/* if all adjacent nodes of j are already visited, pop j from
				 * recursive stack and push j onto output stack */
				head-- ;
				Stack[--top] = j ;
				PRINTF ("  end   dfs at %d ] head : %d\n", j, head) ;
			}
		}

		plength = l_length ;
		return (top) ;
	}

	/**
	 * Finds the pattern of x, for the solution of Lx=b.
	 *
	 * @param n L is n-by-n, where n >= 0
	 * @param k also used as the mark value, for the Flag array
	 * @param Ap
	 * @param Ai
	 * @param Q
	 * @param Pinv Pinv[i] = k if i is kth pivot row, or EMPTY if row i
	 * is not yet pivotal.
	 * @param Stack size n
	 * @param Flag size n.  Initially, all of Flag[0..n-1] < k.  After
	 * lsolve_symbolicis done, Flag[i] == k if i is in
	 * the pattern of the output, and Flag[0..n-1] <= k.
	 * @param Lpend for symmetric pruning
	 * @param Ap_pos workspace used in dfs
	 * @param LU LU factors (pattern and values)
	 * @param lup pointer to free space in LU
	 * @param Llen size n, Llen[k] = # nonzeros in column k of L
	 * @param Lip size n, Lip[k] is position in LU of column k of L
	 * @param k1 the block of A is from k1 to k2-1
	 * @param PSinv inverse of P from symbolic factorization
	 * @return
	 */
	public static int lsolve_symbolic(int n, int k, int[] Ap, int[] Ai,
			int[] Q, int[] Pinv, int[] Stack, int[] Flag, int[] Lpend,
			int[] Ap_pos, double LU, int lup, int[] Llen, int[] Lip,
			int k1, int[] PSinv)
	{
		int[] Lik;
		int i, p, pend, oldcol, kglobal, top, l_length;

		top = n ;
		l_length = 0 ;
		Lik = (int[]) (LU + lup);

		/* ---------------------------------------------------------------------- */
		/* BTF factorization of A (k1:k2-1, k1:k2-1) */
		/* ---------------------------------------------------------------------- */

		kglobal = k + k1 ;  /* column k of the block is col kglobal of A */
		oldcol = Q [kglobal] ;      /* Q must be present for BTF case */
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
			i = PSinv [Ai [p]] - k1 ;
			if (i < 0) continue ;   /* skip entry outside the block */

			/* (i,k) is an entry in the block.  start a DFS at node i */
			PRINTF ("\n ===== DFS at node %d in b, inew: %d\n", i, Pinv [i]) ;
			if (Flag [i] != k)
			{
				if (Pinv [i] >= 0)
				{
					top = dfs (i, k, Pinv, Llen, Lip, Stack, Flag,
							   Lpend, top, LU, Lik, l_length, Ap_pos) ;
				}
				else
				{
					/* i is not pivotal, and not flagged. Flag and put in L */
					Flag [i] = k ;
					Lik [l_length] = i ;
					l_length++;
				}
			}
		}

		/* If Llen [k] is zero, the matrix is structurally singular */
		Llen [k] = l_length ;
		return (top) ;
	}

	/**
	 * Construct the kth column of A, and the off-diagonal part, if requested.
	 * Scatter the numerical values into the workspace X, and construct the
	 * corresponding column of the off-diagonal matrix.
	 *
	 * @param k the column of A (or the column of the block) to get
	 * @param Ap
	 * @param Ai
	 * @param Ax
	 * @param Q column pre-ordering
	 * @param X
	 * @param k1 the block of A is from k1 to k2-1
	 * @param PSinv inverse of P from symbolic factorization
	 * @param Rs scale factors for A
	 * @param scale 0: no scaling, nonzero: scale the rows with Rs
	 * @param Offp off-diagonal matrix (modified by this routine)
	 * @param Offi
	 * @param Offx
	 */
	public static void construct_column(int k, int[] Ap, int[] Ai, double[] Ax,
			int[] Q, double[] X, int k1, int[] PSinv, double[] Rs, int scale,
			int[] Offp, int[] Offi, double[] Offx)
	{
		double aik ;
		int i, p, pend, oldcol, kglobal, poff, oldrow ;

		/* ---------------------------------------------------------------------- */
		/* Scale and scatter the column into X. */
		/* ---------------------------------------------------------------------- */

		kglobal = k + k1 ;          /* column k of the block is col kglobal of A */
		poff = Offp [kglobal] ;     /* start of off-diagonal column */
		oldcol = Q [kglobal] ;
		pend = Ap [oldcol+1] ;

		if (scale <= 0)
		{
			/* no scaling */
			for (p = Ap [oldcol] ; p < pend ; p++)
			{
				oldrow = Ai [p] ;
				i = PSinv [oldrow] - k1 ;
				aik = Ax [p] ;
				if (i < 0)
				{
					/* this is an entry in the off-diagonal part */
					Offi [poff] = oldrow ;
					Offx [poff] = aik ;
					poff++ ;
				}
				else
				{
					/* (i,k) is an entry in the block.  scatter into X */
					X [i] = aik ;
				}
			}
		}
		else
		{
			/* row scaling */
			for (p = Ap [oldcol] ; p < pend ; p++)
			{
				oldrow = Ai [p] ;
				i = PSinv [oldrow] - k1 ;
				aik = Ax [p] ;
				aik /= Rs [oldrow] ;
				//SCALE_DIV (aik, Rs [oldrow]) ;
				if (i < 0)
				{
					/* this is an entry in the off-diagonal part */
					Offi [poff] = oldrow ;
					Offx [poff] = aik ;
					poff++ ;
				}
				else
				{
					/* (i,k) is an entry in the block.  scatter into X */
					X [i] = aik ;
				}
			}
		}

		Offp [kglobal+1] = poff ;   /* start of the next col of off-diag part */
	}

	/**
	 * Computes the numerical values of x, for the solution of Lx=b.  Note that x
	 * may include explicit zeros if numerical cancelation occurs.  L is assumed
	 * to be unit-diagonal, with possibly unsorted columns (but the first entry in
	 * the column must always be the diagonal entry).
	 *
	 * @param Pinv Pinv[i] = k if i is kth pivot row, or EMPTY if row i
	 * is not yet pivotal.
	 * @param LU LU factors (pattern and values)
	 * @param Stack stack for dfs
	 * @param Lip size n, Lip[k] is position in LU of column k of L
	 * @param top top of stack on input
	 * @param n A is n-by-n
	 * @param Llen size n, Llen[k] = # nonzeros in column k of L
	 * @param X size n, initially zero.  On output,
	 * X[Ui[up1..up-1]] and X[Li[lp1..lp-1]] contains the solution.
	 */
	public static void lsolve_numeric(int[] Pinv, double LU, int[] Stack,
			int[] Lip, int top, int n, int[] Llen, double[] X)
	{
		double xj;
		double[] Lx;
		int[] Li;
		int p, s, j, jnew, len;

		/* solve Lx=b */
		for (s = top ; s < n ; s++)
		{
			/* forward solve with column j of L */
			j = Stack [s] ;
			jnew = Pinv [j] ;
			ASSERT (jnew >= 0) ;
			xj = X [j] ;
			GET_POINTER (LU, Lip, Llen, Li, Lx, jnew, len) ;
			ASSERT (Lip [jnew] <= Lip [jnew+1]) ;
			for (p = 0 ; p < len ; p++)
			{
				X [Li [p]] -= Lx [p] * xj ;
				//MULT_SUB (X [Li [p]], Lx [p], xj) ;
			}
		}
	}

	/**
	 * Find a pivot via partial pivoting, and scale the column of L.
	 *
	 * @param diagrow
	 * @param p_pivrow
	 * @param p_pivot
	 * @param p_abs_pivot
	 * @param tol
	 * @param X
	 * @param LU LU factors (pattern and values)
	 * @param Lip
	 * @param Llen
	 * @param k
	 * @param n
	 * @param Pinv Pinv[i] = k if row i is kth pivot row, or EMPTY if
	 * row i is not yet pivotal.
	 * @param p_firstrow
	 * @param Common
	 * @return
	 */
	public static int lpivot(int diagrow, int p_pivrow, double p_pivot,
			double p_abs_pivot, double tol, double[] X, double[] LU, int[] Lip,
			int[] Llen, int k, int n, int[] Pinv , int p_firstrow,
			KLU_common Common)
	{
		double x, pivot;
		double[] Lx;
		double abs_pivot, xabs;
		int p, i, ppivrow, pdiag, pivrow, last_row_index, firstrow, len;
		int[] Li;

		pivrow = EMPTY ;
		if (Llen [k] == 0)
		{
			/* matrix is structurally singular */
			if (Common.halt_if_singular != 0)
			{
				return (FALSE) ;
			}
			for (firstrow = p_firstrow ; firstrow < n ; firstrow++)
			{
				PRINTF ("check %d\n", firstrow) ;
				if (Pinv [firstrow] < 0)
				{
					/* found the lowest-numbered non-pivotal row.  Pick it. */
					pivrow = firstrow ;
					PRINTF ("Got pivotal row: %d\n", pivrow) ;
					break ;
				}
			}
			ASSERT (pivrow >= 0 && pivrow < n) ;
			pivot = 0.0 ;
			//CLEAR (pivot) ;
			p_pivrow = pivrow ;   // FIXME pass by reference
			p_pivot = pivot ;
			p_abs_pivot = 0 ;
			p_firstrow = firstrow ;
			return (FALSE) ;
		}

		pdiag = EMPTY ;
		ppivrow = EMPTY ;
		abs_pivot = EMPTY ;
		i = Llen [k] - 1 ;
		GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
		last_row_index = Li [i] ;

		/* decrement the length by 1 */
		Llen [k] = i ;
		GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;

		/* look in Li [0 ..Llen [k] - 1 ] for a pivot row */
		for (p = 0 ; p < len ; p++)
		{
			/* gather the entry from X and store in L */
			i = Li [p] ;
			x = X [i] ;
			X [i] = 0.0 ;
			//CLEAR (X [i]) ;

			Lx [p] = x ;
			xabs = Math.abs( x ) ;
			//ABS (xabs, x) ;

			/* find the diagonal */
			if (i == diagrow)
			{
				pdiag = p ;
			}

			/* find the partial-pivoting choice */
			if (xabs > abs_pivot)
			{
				abs_pivot = xabs ;
				ppivrow = p ;
			}
		}

		xabs = Math.abs (X [last_row_index]) ;
		//ABS (xabs, X [last_row_index]) ;
		if (xabs > abs_pivot)
		{
			abs_pivot = xabs ;
			ppivrow = EMPTY ;
		}

		/* compare the diagonal with the largest entry */
		if (last_row_index == diagrow)
		{
			if (xabs >= tol * abs_pivot)
			{
				abs_pivot = xabs ;
				ppivrow = EMPTY ;
			}
		}
		else if (pdiag != EMPTY)
		{
			xabs = Math.abs (Lx [pdiag]) ;
			//ABS (xabs, Lx [pdiag]) ;
			if (xabs >= tol * abs_pivot)
			{
				/* the diagonal is large enough */
				abs_pivot = xabs ;
				ppivrow = pdiag ;
			}
		}

		if (ppivrow != EMPTY)
		{
			pivrow = Li [ppivrow] ;
			pivot  = Lx [ppivrow] ;
			/* overwrite the ppivrow values with last index values */
			Li [ppivrow] = last_row_index ;
			Lx [ppivrow] = X [last_row_index] ;
		}
		else
		{
			pivrow = last_row_index ;
			pivot = X [last_row_index] ;
		}
		X [last_row_index] = 0.0 ;
		//CLEAR (X [last_row_index]) ;

		p_pivrow = pivrow ;
		p_pivot = pivot ;
		p_abs_pivot = abs_pivot ;
		ASSERT (pivrow >= 0 && pivrow < n) ;

		if (IS_ZERO (pivot) && Common.halt_if_singular != 0)
		{
			/* numerically singular case */
			return (FALSE) ;
		}

		/* divide L by the pivot value */
		for (p = 0 ; p < Llen [k] ; p++)
		{
			Lx [p] /= pivot ;
			//DIV (Lx [p], Lx [p], pivot) ;
		}

		return (TRUE) ;
	}

	/**
	 * Prune the columns of L to reduce work in subsequent depth-first searches.
	 *
	 * @param Lpend Lpend[j] marks symmetric pruning point for L(:,j)
	 * @param Pinv Pinv[i] = k if row i is kth pivot row, or EMPTY if
	 * row i is not yet pivotal.
	 * @param k pruneusing column k of U
	 * @param pivrow current pivot row
	 * @param LU LU factors (pattern and values)
	 * @param Uip size n, column pointers for U
	 * @param Lip size n, column pointers for L
	 * @param Ulen size n, column length of U
	 * @param Llen size n, column length of L
	 */
	public static void prune(int[] Lpend, int[] Pinv, int k, int pivrow,
			double LU, int[] Uip, int[] Lip, int[] Ulen, int[] Llen)
	{
		double x;
		double[] Lx, Ux;
		int[] Li, Ui;
		int p, i, j, p2, phead, ptail, llen, ulen;

		/* check to see if any column of L can be pruned */
		GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
		for (p = 0 ; p < ulen ; p++)
		{
			j = Ui [p] ;
			ASSERT (j < k) ;
			PRINTF ("%d is pruned: %d. Lpend[j] %d Lip[j+1] %d\n",
				j, Lpend [j] != EMPTY, Lpend [j], Lip [j+1]) ;
			if (Lpend [j] == EMPTY)
			{
				/* scan column j of L for the pivot row */
				GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
				for (p2 = 0 ; p2 < llen ; p2++)
				{
					if (pivrow == Li [p2])
					{
						/* found it!  This column can be pruned */
						if (!NDEBUG)
						{
							PRINTF ("==== PRUNE: col j %d of L\n", j) ;
							{
								int p3 ;
								for (p3 = 0 ; p3 < Llen [j] ; p3++)
								{
									PRINTF ("before: %i  pivotal: %d\n", Li [p3],
												Pinv [Li [p3]] >= 0) ;
								}
							}
						}

						/* partition column j of L.  The unit diagonal of L
						 * is not stored in the column of L. */
						phead = 0 ;
						ptail = Llen [j] ;
						while (phead < ptail)
						{
							i = Li [phead] ;
							if (Pinv [i] >= 0)
							{
								/* leave at the head */
								phead++ ;
							}
							else
							{
								/* swap with the tail */
								ptail-- ;
								Li [phead] = Li [ptail] ;
								Li [ptail] = i ;
								x = Lx [phead] ;
								Lx [phead] = Lx [ptail] ;
								Lx [ptail] = x ;
							}
						}

						/* set Lpend to one past the last entry in the
						 * first part of the column of L.  Entries in
						 * Li [0 ... Lpend [j]-1] are the only part of
						 * column j of L that needs to be scanned in the DFS.
						 * Lpend [j] was EMPTY; setting it >= 0 also flags
						 * column j as pruned. */
						Lpend [j] = ptail ;

						if (!NDEBUG)
						{
							int p3 ;
							for (p3 = 0 ; p3 < Llen [j] ; p3++)
							{
								if (p3 == Lpend [j]) PRINTF (("----\n")) ;
								PRINTF ("after: %i  pivotal: %d\n", Li [p3],
											Pinv [Li [p3]] >= 0) ;
							}
						}

						break ;
					}
				}
			}
		}
	}

	/**
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers for A
	 * @param Ai size nz = Ap[n], row indices for A
	 * @param Ax size nz, values of A
	 * @param Q size n, optional input permutation
	 * @param lusize initial size of LU on input
	 * @param Pinv size n, inverse row permutation, where Pinv[i] = k if
	 * row i is the kth pivot row
	 * @param P size n, row permutation, where P[k] = i if row i is the
	 * kth pivot row.
	 * @param p_LU LU array, size lusize on input
	 * @param Udiag size n, diagonal of U
	 * @param Llen size n, column length of L
	 * @param Ulen size n, column length of U
	 * @param Lip size n, column pointers for L
	 * @param Uip size n, column pointers for U
	 * @param lnz size of L
	 * @param unz size of U
	 * @param X size n, undefined on input, zero on output
	 * @param Stack size n
	 * @param Flag size n
	 * @param Ap_pos size n
	 * @param Lpend size n workspace, for pruning only
	 * @param k1 the block of A is from k1 to k2-1
	 * @param PSinv inverse of P from symbolic factorization
	 * @param Rs scale factors for A
	 * @param Offp off-diagonal matrix (modified by this routine)
	 * @param Offi
	 * @param Offx
	 * @param Common
	 * @return final size of LU on output
	 */
	public static int klu_kernel(int n, int[] Ap, int[] Ai, double[] Ax,
			int[] Q, int lusize, int[] Pinv, int[] P, double[][] p_LU,
			double[] Udiag, int[] Llen, int[] Ulen, int[] Lip, int[] Uip,
			int lnz, int unz, double[] X, int[] Stack, int[] Flag,
			int[] Ap_pos, int[] Lpend, int k1, int[] PSinv, double[] Rs,
			int[] Offp, int[] Offi, double[] Offx, KLU_common Common)
	{
		double pivot;
		double abs_pivot, xsize, nunits, tol, memgrow;
		double[] Ux;
		int[] Li, Ui;
		double[] LU;          /* LU factors (pattern and values) */
		int k, p, i, j, pivrow = 0, kbar, diagrow, firstrow, lup, top, scale, len;
		int newlusize;

		double[] Lx;  // only used when debugging

		ASSERT (Common != null) ;
		scale = Common.scale ;
		tol = Common.tol ;
		memgrow = Common.memgrow ;
		lnz = 0 ;  // FIXME pass by reference
		unz = 0 ;
		pivot = 0.0 ;
		//CLEAR (pivot) ;

		/* ---------------------------------------------------------------------- */
		/* get initial Li, Lx, Ui, and Ux */
		/* ---------------------------------------------------------------------- */

		PRINTF ("input: lusize %d \n", lusize) ;
		ASSERT (lusize > 0) ;
		LU = p_LU ;

		/* ---------------------------------------------------------------------- */
		/* initializations */
		/* ---------------------------------------------------------------------- */

		firstrow = 0 ;
		lup = 0 ;

		for (k = 0 ; k < n ; k++)
		{
			X [k] = 0 ;
			//CLEAR (X [k]) ;
			Flag [k] = EMPTY ;
			Lpend [k] = EMPTY ;     /* flag k as not pruned */
		}

		/* ---------------------------------------------------------------------- */
		/* mark all rows as non-pivotal and determine initial diagonal mapping */
		/* ---------------------------------------------------------------------- */

		/* PSinv does the symmetric permutation, so don't do it here */
		for (k = 0 ; k < n ; k++)
		{
			P [k] = k ;
			Pinv [k] = FLIP (k) ;   /* mark all rows as non-pivotal */
		}
		/* initialize the construction of the off-diagonal matrix */
		Offp [0] = 0 ;

		/* P [k] = row means that UNFLIP (Pinv [row]) = k, and visa versa.
		 * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
		 * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
		 * pivotal. */

		if (!NDEBUG)
		{
			for (k = 0 ; k < n ; k++)
			{
				PRINTF ("Initial P [%d] = %d\n", k, P [k]) ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* factorize */
		/* ---------------------------------------------------------------------- */

		for (k = 0 ; k < n ; k++)
		{

			PRINTF ("\n\n==================================== k: %d\n", k) ;

			/* ------------------------------------------------------------------ */
			/* determine if LU factors have grown too big */
			/* ------------------------------------------------------------------ */

			/* (n - k) entries for L and k entries for U */
			nunits = DUNITS (Integer.class, n - k) + DUNITS (Integer.class, k) +
					DUNITS (Double.class, n - k) + DUNITS (Double.class, k) ;

			/* LU can grow by at most 'nunits' entries if the column is dense */
			PRINTF ("lup %d lusize %g lup+nunits: %g\n", lup, (double) lusize,
				lup+nunits) ;
			xsize = ((double) lup) + nunits ;
			if (xsize > (double) lusize)
			{
				/* check here how much to grow */
				xsize = (memgrow * ((double) lusize) + 4*n + 1) ;
				if (INT_OVERFLOW (xsize))
				{
					PRINTF ("Matrix is too large (int overflow)\n") ;
					Common.status = KLU_TOO_LARGE ;
					return (lusize) ;
				}
				newlusize = (int) (memgrow * lusize + 2*n + 1) ;
				/* Future work: retry mechanism in case of malloc failure */
				LU = klu_realloc (newlusize, lusize, Double.class, LU, Common) ;
				Common.nrealloc++ ;
				p_LU = LU ;
				if (Common.status == KLU_OUT_OF_MEMORY)
				{
					PRINTF ("Matrix is too large (LU)\n") ;
					return (lusize) ;
				}
				lusize = newlusize ;
				PRINTF ("inc LU to %d done\n", lusize) ;
			}

			/* ------------------------------------------------------------------ */
			/* start the kth column of L and U */
			/* ------------------------------------------------------------------ */

			Lip [k] = lup ;

			/* ------------------------------------------------------------------ */
			/* compute the nonzero pattern of the kth column of L and U */
			/* ------------------------------------------------------------------ */

			if (!NDEBUG)
			{
				for (i = 0 ; i < n ; i++)
				{
					ASSERT (Flag [i] < k) ;
					/* ASSERT (X [i] == 0) ; */
					ASSERT (IS_ZERO (X [i])) ;
				}
			}

			top = lsolve_symbolic (n, k, Ap, Ai, Q, Pinv, Stack, Flag,
					Lpend, Ap_pos, LU, lup, Llen, Lip, k1, PSinv) ;

			if (!NDEBUG)
			{
				PRINTF ("--- in U:\n") ;
				for (p = top ; p < n ; p++)
				{
					PRINTF ("pattern of X for U: %d : %d pivot row: %d\n",
						p, Stack [p], Pinv [Stack [p]]) ;
					ASSERT (Flag [Stack [p]] == k) ;
				}
				PRINTF ("--- in L:\n") ;
				Li = (int[]) (LU + Lip [k]);
				for (p = 0 ; p < Llen [k] ; p++)
				{
					PRINTF ("pattern of X in L: %d : %d pivot row: %d\n",
						p, Li [p], Pinv [Li [p]]) ;
					ASSERT (Flag [Li [p]] == k) ;
				}
				p = 0 ;
				for (i = 0 ; i < n ; i++)
				{
					ASSERT (Flag [i] <= k) ;
					if (Flag [i] == k) p++ ;
				}
			}

			/* ------------------------------------------------------------------ */
			/* get the column of the matrix to factorize and scatter into X */
			/* ------------------------------------------------------------------ */

			construct_column (k, Ap, Ai, Ax, Q, X,
				k1, PSinv, Rs, scale, Offp, Offi, Offx) ;

			/* ------------------------------------------------------------------ */
			/* compute the numerical values of the kth column (s = L \ A (:,k)) */
			/* ------------------------------------------------------------------ */

			lsolve_numeric (Pinv, LU, Stack, Lip, top, n, Llen, X) ;

			if (!NDEBUG)
			{
				for (p = top ; p < n ; p++)
				{
					PRINTF ("X for U %d : ",  Stack [p]) ;
					PRINT_ENTRY (X [Stack [p]]) ;
				}
				Li = (int[]) (LU + Lip [k]) ;
				for (p = 0 ; p < Llen [k] ; p++)
				{
					PRINTF ("X for L %d : ", Li [p]) ;
					PRINT_ENTRY (X [Li [p]]) ;
				}
			}

			/* ------------------------------------------------------------------ */
			/* partial pivoting with diagonal preference */
			/* ------------------------------------------------------------------ */

			/* determine what the "diagonal" is */
			diagrow = P [k] ;   /* might already be pivotal */
			PRINTF ("k %d, diagrow = %d, UNFLIP (diagrow) = %d\n",
				k, diagrow, UNFLIP (diagrow)) ;

			/* find a pivot and scale the pivot column */
			if (lpivot (diagrow, pivrow, pivot, abs_pivot, tol, X, LU, Lip,
						Llen, k, n, Pinv, firstrow, Common) == 0)
			{
				/* matrix is structurally or numerically singular */
				Common.status = KLU_SINGULAR ;
				if (Common.numerical_rank == EMPTY)
				{
					Common.numerical_rank = k+k1 ;
					Common.singular_col = Q [k+k1] ;
				}
				if (Common.halt_if_singular != 0)
				{
					/* do not continue the factorization */
					return (lusize) ;
				}
			}

			/* we now have a valid pivot row, even if the column has NaN's or
			 * has no entries on or below the diagonal at all. */
			PRINTF ("\nk %d : Pivot row %d : ", k, pivrow) ;
			PRINT_ENTRY (pivot) ;
			ASSERT (pivrow >= 0 && pivrow < n) ;
			ASSERT (Pinv [pivrow] < 0) ;

			/* set the Uip pointer */
			Uip [k] = Lip [k] + UNITS (Integer.class, Llen [k]) + UNITS (Double.class, Llen [k]) ;

			/* move the lup pointer to the position where indices of U
			 * should be stored */
			lup += UNITS (Integer.class, Llen [k]) + UNITS (Double.class, Llen [k]) ;

			Ulen [k] = n - top ;

			/* extract Stack [top..n-1] to Ui and the values to Ux and clear X */
			GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
			for (p = top, i = 0 ; p < n ; p++, i++)
			{
				j = Stack [p] ;
				Ui [i] = Pinv [j] ;
				Ux [i] = X [j] ;
				X [j] = 0.0 ;
				//CLEAR (X [j]) ;
			}

			/* position the lu index at the starting point for next column */
			lup += UNITS (Integer.class, Ulen [k]) + UNITS (Double.class, Ulen [k]) ;

			/* U(k,k) = pivot */
			Udiag [k] = pivot ;

			/* ------------------------------------------------------------------ */
			/* log the pivot permutation */
			/* ------------------------------------------------------------------ */

			ASSERT (UNFLIP (Pinv [diagrow]) < n) ;
			ASSERT (P [UNFLIP (Pinv [diagrow])] == diagrow) ;

			if (pivrow != diagrow)
			{
				/* an off-diagonal pivot has been chosen */
				Common.noffdiag++ ;
				PRINTF (">>>>>>>>>>>>>>>>> pivrow %d k %d off-diagonal\n",
							pivrow, k) ;
				if (Pinv [diagrow] < 0)
				{
					/* the former diagonal row index, diagrow, has not yet been
					 * chosen as a pivot row.  Log this diagrow as the "diagonal"
					 * entry in the column kbar for which the chosen pivot row,
					 * pivrow, was originally logged as the "diagonal" */
					kbar = FLIP (Pinv [pivrow]) ;
					P [kbar] = diagrow ;
					Pinv [diagrow] = FLIP (kbar) ;
				}
			}
			P [k] = pivrow ;
			Pinv [pivrow] = k ;

			if (!NDEBUG)
			{
				for (i = 0 ; i < n ; i++) { ASSERT (IS_ZERO (X [i])) ;}
				GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
				for (p = 0 ; p < len ; p++)
				{
					PRINTF ("Column %d of U: %d : ", k, Ui [p]) ;
					PRINT_ENTRY (Ux [p]) ;
				}
				GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
				for (p = 0 ; p < len ; p++)
				{
					PRINTF ("Column %d of L: %d : ", k, Li [p]) ;
					PRINT_ENTRY (Lx [p]) ;
				}
			}

			/* ------------------------------------------------------------------ */
			/* symmetric pruning */
			/* ------------------------------------------------------------------ */

			prune (Lpend, Pinv, k, pivrow, LU, Uip, Lip, Ulen, Llen) ;

			lnz += Llen [k] + 1 ; /* 1 added to lnz for diagonal */
			unz += Ulen [k] + 1 ; /* 1 added to unz for diagonal */
		}

		/* ---------------------------------------------------------------------- */
		/* finalize column pointers for L and U, and put L in the pivotal order */
		/* ---------------------------------------------------------------------- */

		for (p = 0 ; p < n ; p++)
		{
			Li = (int[]) (LU + Lip [p]) ;
			for (i = 0 ; i < Llen [p] ; i++)
			{
				Li [i] = Pinv [Li [i]] ;
			}
		}

		if (!NDEBUG)
		{
			for (i = 0 ; i < n ; i++)
			{
				PRINTF ("P [%d] = %d   Pinv [%d] = %d\n", i, P [i], i, Pinv [i]) ;
			}
			for (i = 0 ; i < n ; i++)
			{
				ASSERT (Pinv [i] >= 0 && Pinv [i] < n) ;
				ASSERT (P [i] >= 0 && P [i] < n) ;
				ASSERT (P [Pinv [i]] == i) ;
				ASSERT (IS_ZERO (X [i])) ;
			}
		}

		/* ---------------------------------------------------------------------- */
		/* shrink the LU factors to just the required size */
		/* ---------------------------------------------------------------------- */

		newlusize = lup ;
		ASSERT ((int) newlusize <= lusize) ;

		/* this cannot fail, since the block is descreasing in size */
		LU = klu_realloc (newlusize, lusize, Double.class, LU, Common) ;
		p_LU = LU ;
		return (newlusize) ;
	}

}
