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
import edu.ufl.cise.klu.common.KLU_symbolic;

/**
 * Given an input permutation P and Q, create the Symbolic object.  BTF can
 * be done to modify the user's P and Q (does not perform the max transversal;
 * just finds the strongly-connected components).
 */
public class Dklu_analyze_given extends Dklu_internal
{

	/**
	 * Allocate Symbolic object, and check input matrix.  Not user callable.
	 *
	 * @param n
	 * @param Ap
	 * @param Ai
	 * @param Common
	 * @return
	 */
	protected static KLU_symbolic klu_alloc_symbolic(int n, int[] Ap, int[] Ai,
		    KLU_common Common)
	{
	    KLU_symbolic Symbolic ;
	    int[] P, Q, R ;
	    double[] Lnz ;
	    int nz, i, j, p, pend ;

	    if (Common == null)
	    {
	        return (null) ;
	    }
	    Common.status = KLU_common.KLU_OK ;

	    /* A is n-by-n, with n > 0.  Ap [0] = 0 and nz = Ap [n] >= 0 required.
	     * Ap [j] <= Ap [j+1] must hold for all j = 0 to n-1.  Row indices in Ai
	     * must be in the range 0 to n-1, and no duplicate entries can be present.
	     * The list of row indices in each column of A need not be sorted.
	     */

	    if (n <= 0 || Ap == null || Ai == null)
	    {
	        /* Ap and Ai must be present, and n must be > 0 */
	        Common.status = KLU_common.KLU_INVALID ;
	        return (null) ;
	    }

	    nz = Ap [n] ;
	    if (Ap [0] != 0 || nz < 0)
	    {
	        /* nz must be >= 0 and Ap [0] must equal zero */
	        Common.status = KLU_common.KLU_INVALID ;
	        return (null) ;
	    }

	    for (j = 0 ; j < n ; j++)
	    {
	        if (Ap [j] > Ap [j+1])
	        {
	            /* column pointers must be non-decreasing */
	            Common.status = KLU_common.KLU_INVALID ;
	            return (null) ;
	        }
	    }
	    P = Dklu_memory.klu_malloc (n, sizeof (Int), Common) ;
	    if (Common.status < KLU_common.KLU_OK)
	    {
	        /* out of memory */
	        Common.status = KLU_common.KLU_OUT_OF_MEMORY ;
	        return (null) ;
	    }
	    for (i = 0 ; i < n ; i++)
	    {
	        P [i] = EMPTY ;
	    }
	    for (j = 0 ; j < n ; j++)
	    {
	        pend = Ap [j+1] ;
	        for (p = Ap [j] ; p < pend ; p++)
	        {
	            i = Ai [p] ;
	            if (i < 0 || i >= n || P [i] == j)
	            {
	                /* row index out of range, or duplicate entry */
	                Dklu_memory.klu_free (P, n, sizeof (Int), Common) ;
	                Common.status = KLU_common.KLU_INVALID ;
	                return (null) ;
	            }
	            /* flag row i as appearing in column j */
	            P [i] = j ;
	        }
	    }

	    /* ---------------------------------------------------------------------- */
	    /* allocate the Symbolic object */
	    /* ---------------------------------------------------------------------- */

	    Symbolic = Dklu_memory.klu_malloc (sizeof (KLU_symbolic), 1, Common) ;
	    if (Common.status < KLU_common.KLU_OK)
	    {
	        /* out of memory */
	        Dklu_memory.klu_free (P, n, sizeof (Int), Common) ;
	        Common.status = KLU_common.KLU_OUT_OF_MEMORY ;
	        return (null) ;
	    }

	    Q = Dklu_memory.klu_malloc (n, sizeof (Int), Common) ;
	    R = Dklu_memory.klu_malloc (n+1, sizeof (Int), Common) ;
	    Lnz = Dklu_memory.klu_malloc (n, sizeof (double), Common) ;

	    Symbolic.n = n ;
	    Symbolic.nz = nz ;
	    Symbolic.P = P ;
	    Symbolic.Q = Q ;
	    Symbolic.R = R ;
	    Symbolic.Lnz = Lnz ;

	    if (Common.status < KLU_common.KLU_OK)
	    {
	        /* out of memory */
	        Dklu_free_symbolic.klu_free_symbolic (Symbolic, Common) ;
	        Common.status = KLU_common.KLU_OUT_OF_MEMORY ;
	        return (null) ;
	    }

	    return (Symbolic) ;
	}

	/**
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Puser size n, user's row permutation (may be null)
	 * @param Quser size n, user's column permutation (may be null)
	 * @param Common
	 * @return
	 */
	public static KLU_symbolic klu_analyze_given(int n, int[] Ap, int[] Ai,
		    int[] Puser, int[] Quser, KLU_common Common)
	{
		KLU_symbolic Symbolic ;
	    double[] Lnz ;
	    int nblocks, nz, block, maxblock, nzoff, p, pend, do_btf, k ;
	    int[] P, Q, R ;

	    /* ---------------------------------------------------------------------- */
	    /* determine if input matrix is valid, and get # of nonzeros */
	    /* ---------------------------------------------------------------------- */

	    Symbolic = Dklu_alloc_symbolic.klu_alloc_symbolic (n, Ap, Ai, Common) ;
	    if (Symbolic == null)
	    {
	        return (null) ;
	    }
	    P = Symbolic.P ;
	    Q = Symbolic.Q ;
	    R = Symbolic.R ;
	    Lnz = Symbolic.Lnz ;
	    nz = Symbolic.nz ;

	    /* ---------------------------------------------------------------------- */
	    /* Q = Quser, or identity if Quser is null */
	    /* ---------------------------------------------------------------------- */

	    if (Quser == null)
	    {
	        for (k = 0 ; k < n ; k++)
	        {
	            Q [k] = k ;
	        }
	    }
	    else
	    {
	        for (k = 0 ; k < n ; k++)
	        {
	            Q [k] = Quser [k] ;
	        }
	    }

	    /* ---------------------------------------------------------------------- */
	    /* get the control parameters for BTF and ordering method */
	    /* ---------------------------------------------------------------------- */

	    do_btf = Common.btf ;
	    do_btf = (do_btf) ? TRUE : FALSE ;
	    Symbolic.ordering = 2 ;
	    Symbolic.do_btf = do_btf ;

	    /* ---------------------------------------------------------------------- */
	    /* find the block triangular form, if requested */
	    /* ---------------------------------------------------------------------- */

	    if (do_btf == 1)
	    {

	        /* ------------------------------------------------------------------ */
	        /* get workspace for BTF_strongcomp */
	        /* ------------------------------------------------------------------ */

	        int[] Pinv, Work, Bi ;
	        int k1, k2, nk, oldcol ;

	        Work = Dklu_memory.klu_malloc (4*n, sizeof (Int), Common) ;
	        Pinv = Dklu_memory.klu_malloc (n, sizeof (Int), Common) ;
	        if (Puser != null)
	        {
	            Bi = Dklu_memory.klu_malloc (nz+1, sizeof (Int), Common) ;
	        }
	        else
	        {
	            Bi = Ai ;
	        }

	        if (Common.status < KLU_common.KLU_OK)
	        {
	            /* out of memory */
	            Dklu_memory.klu_free (Work, 4*n, sizeof (Int), Common) ;
	            Dklu_memory.klu_free (Pinv, n, sizeof (Int), Common) ;
	            if (Puser != null)
	            {
	                Dklu_memory.klu_free (Bi, nz+1, sizeof (Int), Common) ;
	            }
	            Dklu_free_symbolic.klu_free_symbolic (Symbolic, Common) ;
	            Common.status = KLU_common.KLU_OUT_OF_MEMORY ;
	            return (null) ;
	        }

	        /* ------------------------------------------------------------------ */
	        /* B = Puser * A */
	        /* ------------------------------------------------------------------ */

	        if (Puser != null)
	        {
	            for (k = 0 ; k < n ; k++)
	            {
	                Pinv [Puser [k]] = k ;
	            }
	            for (p = 0 ; p < nz ; p++)
	            {
	                Bi [p] = Pinv [Ai [p]] ;
	            }
	        }

	        /* ------------------------------------------------------------------ */
	        /* find the strongly-connected components */
	        /* ------------------------------------------------------------------ */

	        /* modifies Q, and determines P and R */
	        nblocks = Dbtf_strongcomp.btf_strongcomp (n, Ap, Bi, Q, P, R, Work) ;

	        /* ------------------------------------------------------------------ */
	        /* P = P * Puser */
	        /* ------------------------------------------------------------------ */

	        if (Puser != null)
	        {
	            for (k = 0 ; k < n ; k++)
	            {
	                Work [k] = Puser [P [k]] ;
	            }
	            for (k = 0 ; k < n ; k++)
	            {
	                P [k] = Work [k] ;
	            }
	        }

	        /* ------------------------------------------------------------------ */
	        /* Pinv = inverse of P */
	        /* ------------------------------------------------------------------ */

	        for (k = 0 ; k < n ; k++)
	        {
	            Pinv [P [k]] = k ;
	        }

	        /* ------------------------------------------------------------------ */
	        /* analyze each block */
	        /* ------------------------------------------------------------------ */

	        nzoff = 0 ;         /* nz in off-diagonal part */
	        maxblock = 1 ;      /* size of the largest block */

	        for (block = 0 ; block < nblocks ; block++)
	        {

	            /* -------------------------------------------------------------- */
	            /* the block is from rows/columns k1 to k2-1 */
	            /* -------------------------------------------------------------- */

	            k1 = R [block] ;
	            k2 = R [block+1] ;
	            nk = k2 - k1 ;
	            PRINTF (("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk)) ;
	            maxblock = MAX (maxblock, nk) ;

	            /* -------------------------------------------------------------- */
	            /* scan the kth block, C */
	            /* -------------------------------------------------------------- */

	            for (k = k1 ; k < k2 ; k++)
	            {
	                oldcol = Q [k] ;
	                pend = Ap [oldcol+1] ;
	                for (p = Ap [oldcol] ; p < pend ; p++)
	                {
	                    if (Pinv [Ai [p]] < k1)
	                    {
	                        nzoff++ ;
	                    }
	                }
	            }

	            /* fill-in not estimated */
	            Lnz [block] = EMPTY ;
	        }

	        /* ------------------------------------------------------------------ */
	        /* free all workspace */
	        /* ------------------------------------------------------------------ */

	        Dklu_memory.klu_free (Work, 4*n, sizeof (Int), Common) ;
	        Dklu_memory.klu_free (Pinv, n, sizeof (Int), Common) ;
	        if (Puser != null)
	        {
	            Dklu_memory.klu_free (Bi, nz+1, sizeof (Int), Common) ;
	        }

	    }
	    else
	    {

	        /* ------------------------------------------------------------------ */
	        /* BTF not requested */
	        /* ------------------------------------------------------------------ */

	        nzoff = 0 ;
	        nblocks = 1 ;
	        maxblock = n ;
	        R [0] = 0 ;
	        R [1] = n ;
	        Lnz [0] = EMPTY ;

	        /* ------------------------------------------------------------------ */
	        /* P = Puser, or identity if Puser is null */
	        /* ------------------------------------------------------------------ */

	        for (k = 0 ; k < n ; k++)
	        {
	            P [k] = (Puser == null) ? k : Puser [k] ;
	        }
	    }

	    /* ---------------------------------------------------------------------- */
	    /* return the symbolic object */
	    /* ---------------------------------------------------------------------- */

	    Symbolic.nblocks = nblocks ;
	    Symbolic.maxblock = maxblock ;
	    Symbolic.lnz = EMPTY ;
	    Symbolic.unz = EMPTY ;
	    Symbolic.nzoff = nzoff ;

	    return (Symbolic) ;
	}

}
