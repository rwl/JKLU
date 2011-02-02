/* KLU: a sparse LU factorization algorithm.
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

package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.tdcomplex.util.Entry;
import edu.ufl.cise.klu.tdcomplex.util.Unit;
import edu.ufl.cise.klu.tdcomplex.util.Util;

/**
 * KLU: factorizes P*A into L*U, using the Gilbert-Peierls algorithm [1], with
 * optional symmetric pruning by Eisenstat and Liu [2].  The code is by Tim
 * Davis.  This algorithm is what appears as the default sparse LU routine in
 * MATLAB version 6.0, and still appears in MATLAB 6.5 as [L,U,P] = lu (A).
 * Note that no column ordering is provided (see COLAMD or AMD for suitable
 * orderings).  SuperLU is based on this algorithm, except that it adds the
 * use of dense matrix operations on "supernodes" (adjacent columns with
 * identical).  This code doesn't use supernodes, thus its name ("Kent" LU,
 * as in "Clark Kent", in contrast with Super-LU...).  This algorithm is slower
 * than SuperLU and UMFPACK for large matrices with lots of nonzeros in their
 * factors (such as for most finite-element problems).  However, for matrices
 * with very sparse LU factors, this algorithm is typically faster than both
 * SuperLU and UMFPACK, since in this case there is little chance to exploit
 * dense matrix kernels (the BLAS).
 *
 * Only one block of A is factorized, in the BTF form.  The input n is the
 * size of the block; k1 is the first row and column in the block.
 *
 * NOTE: no error checking is done on the inputs.  This version is not meant to
 * be called directly by the user.  Use klu_factor instead.
 *
 * No fill-reducing ordering is provided.  The ordering quality of
 * klu_kernel_factor is the responsibility of the caller.  The input A must
 * pre-permuted to reduce fill-in, or fill-reducing input permutation Q must
 * be provided.
 *
 * The input matrix A must be in compressed-column form, with either sorted
 * or unsorted row indices.  Row indices for column j of A is in
 * Ai [Ap [j] ... Ap [j+1]-1] and the same range of indices in Ax holds the
 * numerical values.  No duplicate entries are allowed.
 *
 * Copyright 2004-2009, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers [1,2].
 *
 * If an input permutation Q is provided, the factorization L*U = A (P,Q)
 * is computed, where P is determined by partial pivoting, and Q is the input
 * ordering.  If the pivot tolerance is less than 1, the "diagonal" entry that
 * KLU attempts to choose is the diagonal of A (Q,Q).  In other words, the
 * input permutation is applied symmetrically to the input matrix.  The output
 * permutation P includes both the partial pivoting ordering and the input
 * permutation.  If Q is NULL, then it is assumed to be the identity
 * permutation.  Q is not modified.
 *
 * [1] Gilbert, J. R. and Peierls, T., "Sparse Partial Pivoting in Time
 *      Proportional to Arithmetic Operations," SIAM J. Sci. Stat. Comp.,
 *      vol 9, pp.  862-874, 1988.
 * [2] Eisenstat, S. C. and Liu, J. W. H., "Exploiting Structural Symmetry in
 *      Unsymmetric Sparse Symbolic Factorization," SIAM J. Matrix Analysis &
 *      Applic., vol 13, pp.  202-211, 1992.
 */
public class DZklu extends Util {
	
	/* @return: 0 if failure, size of LU if OK */
	public static int KLU_kernel_factor(
			/* inputs, not modified */
			int n,          /* A is n-by-n. n must be > 0. */
			int [ ] Ap,     /* size n+1, column pointers for A */
			int [ ] Ai,     /* size nz = Ap [n], row indices for A */
			Entry [ ] Ax,   /* size nz, values of A */
			int [ ] Q,      /* size n, optional column permutation */
			double Lsize,   /* estimate of number of nonzeros in L */

			/* outputs, not defined on input */
			Unit p_LU,          /* row indices and values of L and U */
			Entry [ ] Udiag,    /* size n, diagonal of U */
			int [ ] Llen,       /* size n, column length of L */
			int [ ] Ulen,       /* size n, column length of U */
			int [ ] Lip,        /* size n, column pointers for L */
			int [ ] Uip,        /* size n, column pointers for U */
			int [ ] P,          /* row permutation, size n */
			int lnz,            /* size of L */
			int unz,            /* size of U */

			/* workspace, undefined on input */
			Entry X,        /* size n double's, zero on output */
			int Work,       /* size 5n int's */

			/* inputs, not modified on output */
			int k1,             /* the block of A is from k1 to k2-1 */
			int [ ] PSinv,      /* inverse of P from symbolic factorization */
			double [ ] Rs,      /* scale factors for A */

			/* inputs, modified on output */
			int [ ] Offp,   /* off-diagonal matrix (modified by this routine) */
			int [ ] Offi,
			Entry [ ] Offx,
			/* --------------- */
			KLU_common Common 
	) {
		double maxlnz, dunits ;
		Unit LU ;
		int Pinv, Lpend, Stack, Flag, Ap_pos, W ;
		int lsize, usize, anz;//, ok ;
		boolean ok ;
//		NativeLong lusize ;
		int lusize ;
		assert (Common != null) ;

		/* ---------------------------------------------------------------------- */
		/* get control parameters, or use defaults */
		/* ---------------------------------------------------------------------- */

		n = MAX (1, n) ;
		anz = Ap [n+k1] - Ap [k1] ;

		if (Lsize <= 0)
		{
			Lsize = -Lsize ;
			Lsize = MAX (Lsize, 1.0) ;
			lsize = (int) (Lsize * anz + n) ;
		}
		else
		{
			lsize = (int) Lsize ;
		}

		usize = lsize ;

		lsize  = MAX (n+1, lsize) ;
		usize  = MAX (n+1, usize) ;

		maxlnz = (((double) n) * ((double) n) + ((double) n)) / 2. ;
		maxlnz = MIN (maxlnz, ((double) INT_MAX)) ;
		lsize  = (int) MIN (maxlnz, lsize) ;
		usize  = (int) MIN (maxlnz, usize) ;

		PRINTF ("Welcome to klu: n %d anz %d k1 %d lsize %d usize %d maxlnz %g\n",
			n, anz, k1, lsize, usize, maxlnz) ;

		/* ---------------------------------------------------------------------- */
		/* allocate workspace and outputs */
		/* ---------------------------------------------------------------------- */

		/* return arguments are not yet assigned */
		p_LU = null ;

		/* these computations are safe from size_t overflow */
		W = Work ;
		Pinv = W ;      W += n ;
		Stack = W ;     W += n ;
		Flag = W ;      W += n ;
		Lpend = W ;     W += n ;
		Ap_pos = W ;    W += n ;

//		dunits = DUNITS (Int, lsize) + DUNITS (Entry, lsize) +
//				DUNITS (Int, usize) + DUNITS (Entry, usize) ;
//		lusize = (size_t) dunits ;
//		ok = !INT_OVERFLOW (dunits) ; 
//		LU = ok ? KLU_malloc (lusize, sizeof (Unit), Common) : null ;
//		if (LU == null)
//		{
//			/* out of memory, or problem too large */
//			Common.status = KLU_OUT_OF_MEMORY ; 
//			lusize = 0 ;
//			return (lusize) ;
//		}
		LU = new Unit();//KLU_malloc (lusize, Unit.class, Common);

		/* ---------------------------------------------------------------------- */
		/* factorize */
		/* ---------------------------------------------------------------------- */

		/* with pruning, and non-recursive depth-first-search */
		lusize = KLU_kernel (n, Ap, Ai, Ax, Q, lusize,
				Pinv, P, LU, Udiag, Llen, Ulen, Lip, Uip, lnz, unz,
				X, Stack, Flag, Ap_pos, Lpend,
				k1, PSinv, Rs, Offp, Offi, Offx, Common) ;

		/* ---------------------------------------------------------------------- */
		/* return LU factors, or return nothing if an error occurred */
		/* ---------------------------------------------------------------------- */

		if (Common.status < KLU_OK)
		{
			LU = null;//KLU_free (LU, lusize, sizeof (Unit), Common) ;
			lusize = 0 ;
		}
		p_LU = LU ;
		PRINTF (" in klu noffdiag %d\n", Common.noffdiag) ;
		return (lusize) ;
	}

}
