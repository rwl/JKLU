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
 * Linear algebraic diagnostics.
 */
public class Dklu_diagnostics extends Dklu_internal
{

	/**
	 * Compute the reciprocal pivot growth factor.  In MATLAB notation:
	 *
	 *   rgrowth = min (max (abs ((R \ A (p,q)) - F))) ./ max (abs (U)))
	 *
	 * Takes O(|A|+|U|) time.
	 *
	 * @param Ap
	 * @param Ai
	 * @param Ax
	 * @param Symbolic
	 * @param Numeric
	 * @param Common
	 * @return TRUE if successful, FALSE otherwise
	 */
	public static int klu_rgrowth(int[] Ap, int[] Ai, double[] Ax,
		    KLU_symbolic Symbolic, KLU_numeric Numeric, KLU_common Common)
	{
		double temp, max_ai, max_ui, min_block_rgrowth ;
	    Entry aik ;
	    int[] Q, Ui, Uip, Ulen, Pinv ;
	    Unit[] LU ;
	    Entry[] Aentry, Ux, Ukk ;
	    double[] Rs ;
	    int i, newrow, oldrow, k1, k2, nk, j, oldcol, k, pend, len ;

	    /* ---------------------------------------------------------------------- */
	    /* check inputs */
	    /* ---------------------------------------------------------------------- */

	    if (Common == null)
	    {
	        return (FALSE) ;
	    }

	    if (Symbolic == null || Ap == null || Ai == null || Ax == null)
	    {
	        Common.status = KLU_common.KLU_INVALID ;
	        return (FALSE) ;
	    }

	    if (Numeric == null)
	    {
	        /* treat this as a singular matrix */
	        Common.rgrowth = 0 ;
	        Common.status = KLU_common.KLU_SINGULAR ;
	        return (TRUE) ;
	    }
	    Common.status = KLU_common.KLU_OK ;

	    /* ---------------------------------------------------------------------- */
	    /* compute the reciprocal pivot growth */
	    /* ---------------------------------------------------------------------- */

	    Aentry = (Entry[]) Ax ;
	    Pinv = Numeric.Pinv ;
	    Rs = Numeric.Rs ;
	    Q = Symbolic.Q ;
	    Common.rgrowth = 1 ;

	    for (i = 0 ; i < Symbolic.nblocks ; i++)
	    {
	        k1 = Symbolic.R[i] ;
	        k2 = Symbolic.R[i+1] ;
	        nk = k2 - k1 ;
	        if (nk == 1)
	        {
	            continue ;      /* skip singleton blocks */
	        }
	        LU = (Unit[]) Numeric.LUbx[i] ;
	        Uip = Numeric.Uip + k1 ;
	        Ulen = Numeric.Ulen + k1 ;
	        Ukk = ((Entry[]) Numeric.Udiag) + k1 ;
	        min_block_rgrowth = 1 ;
	        for (j = 0 ; j < nk ; j++)
	        {
	            max_ai = 0 ;
	            max_ui = 0 ;
	            oldcol = Q[j + k1] ;
	            pend = Ap [oldcol + 1] ;
	            for (k = Ap [oldcol] ; k < pend ; k++)
	            {
	                oldrow = Ai [k] ;
	                newrow = Pinv [oldrow] ;
	                if (newrow < k1)
	                {
	                    continue ;  /* skip entry outside the block */
	                }
	                ASSERT (newrow < k2) ;
	                if (Rs != null)
	                {
	                    /* aik = Aentry [k] / Rs [oldrow] */
	                    SCALE_DIV_ASSIGN (aik, Aentry [k], Rs [newrow]) ;
	                }
	                else
	                {
	                    aik = Aentry [k] ;
	                }
	                /* temp = ABS (aik) */
	                ABS (temp, aik) ;
	                if (temp > max_ai)
	                {
	                    max_ai = temp ;
	                }
	            }

	            GET_POINTER (LU, Uip, Ulen, Ui, Ux, j, len) ;
	            for (k = 0 ; k < len ; k++)
	            {
	                /* temp = ABS (Ux [k]) */
	                ABS (temp, Ux [k]) ;
	                if (temp > max_ui)
	                {
	                    max_ui = temp ;
	                }
	            }
	            /* consider the diagonal element */
	            ABS (temp, Ukk [j]) ;
	            if (temp > max_ui)
	            {
	                max_ui = temp ;
	            }

	            /* if max_ui is 0, skip the column */
	            if (SCALAR_IS_ZERO (max_ui))
	            {
	                continue ;
	            }
	            temp = max_ai / max_ui ;
	            if (temp < min_block_rgrowth)
	            {
	                min_block_rgrowth = temp ;
	            }
	        }

	        if (min_block_rgrowth < Common.rgrowth)
	        {
	            Common.rgrowth = min_block_rgrowth ;
	        }
	    }
	    return (TRUE) ;
	}

	/**
	 * Estimate the condition number.  Uses Higham and Tisseur's algorithm
	 * (A block algorithm for matrix 1-norm estimation, with applications to
	 * 1-norm pseudospectra, SIAM J. Matrix Anal. Appl., 21(4):1185-1201, 2000.
	 *
	 * Takes about O(|A|+5*(|L|+|U|)) time
	 *
	 * @param Ap
	 * @param Ax
	 * @param Symbolic
	 * @param Numeric
	 * @param Common
	 * @return TRUE if successful, FALSE otherwise
	 */
	public static int klu_condest(int[] Ap, double[] Ax, KLU_symbolic Symbolic,
		    KLU_numeric Numeric, KLU_common Common)
	{
		double xj, Xmax, csum, anorm, ainv_norm, est_old, est_new, abs_value ;
	    Entry[] Udiag, Aentry, X, S ;
	    int[] R ;
	    int nblocks, i, j, jmax, jnew, pend, n ;
	    int unchanged ;

	    /* ---------------------------------------------------------------------- */
	    /* check inputs */
	    /* ---------------------------------------------------------------------- */

	    if (Common == null)
	    {
	        return (FALSE) ;
	    }
	    if (Symbolic == null || Ap == null || Ax == null)
	    {
	        Common.status = KLU_common.KLU_INVALID ;
	        return (FALSE) ;
	    }
	    abs_value = 0 ;
	    if (Numeric == null)
	    {
	        /* treat this as a singular matrix */
	        Common.condest = 1 / abs_value ;
	        Common.status = KLU_common.KLU_SINGULAR ;
	        return (TRUE) ;
	    }
	    Common.status = KLU_common.KLU_OK ;

	    /* ---------------------------------------------------------------------- */
	    /* get inputs */
	    /* ---------------------------------------------------------------------- */

	    n = Symbolic.n ;
	    nblocks = Symbolic.nblocks ;
	    R = Symbolic.R ;
	    Udiag = Numeric.Udiag ;

	    /* ---------------------------------------------------------------------- */
	    /* check if diagonal of U has a zero on it */
	    /* ---------------------------------------------------------------------- */

	    for (i = 0 ; i < n ; i++)
	    {
	        ABS (abs_value, Udiag [i]) ;
	        if (SCALAR_IS_ZERO (abs_value))
	        {
	            Common.condest = 1 / abs_value ;
	            Common.status = KLU_common.KLU_SINGULAR ;
	            return (TRUE) ;
	        }
	    }

	    /* ---------------------------------------------------------------------- */
	    /* compute 1-norm (maximum column sum) of the matrix */
	    /* ---------------------------------------------------------------------- */

	    anorm =  0.0 ;
	    Aentry = (Entry[]) Ax ;
	    for (i = 0 ; i < n ; i++)
	    {
	        pend = Ap [i + 1] ;
	        csum = 0.0 ;
	        for (j = Ap [i] ; j < pend ; j++)
	        {
	            ABS (abs_value, Aentry [j]) ;
	            csum += abs_value ;
	        }
	        if (csum > anorm)
	        {
	            anorm = csum ;
	        }
	    }

	    /* ---------------------------------------------------------------------- */
	    /* compute estimate of 1-norm of inv (A) */
	    /* ---------------------------------------------------------------------- */

	    /* get workspace (size 2*n Entry's) */
	    X = Numeric.Xwork ;            /* size n space used in KLU_solve, tsolve */
	    X += n ;                        /* X is size n */
	    S = X + n ;                     /* S is size n */

	    for (i = 0 ; i < n ; i++)
	    {
	        CLEAR (S [i]) ;
	        CLEAR (X [i]) ;
	        REAL (X [i]) = 1.0 / ((double) n) ;
	    }
	    jmax = 0 ;

	    ainv_norm = 0.0 ;
	    for (i = 0 ; i < 5 ; i++)
	    {
	        if (i > 0)
	        {
	            /* X [jmax] is the largest entry in X */
	            for (j = 0 ; j < n ; j++)
	            {
	                /* X [j] = 0 ;*/
	                CLEAR (X [j]) ;
	            }
	            REAL (X [jmax]) = 1 ;
	        }

	        Dklu_solve.klu_solve (Symbolic, Numeric, n, 1, (double[]) X, Common) ;
	        est_old = ainv_norm ;
	        ainv_norm = 0.0 ;

	        for (j = 0 ; j < n ; j++)
	        {
	            /* ainv_norm += ABS (X [j]) ;*/
	            ABS (abs_value, X [j]) ;
	            ainv_norm += abs_value ;
	        }

	        unchanged = TRUE ;

	        for (j = 0 ; j < n ; j++)
	        {
	            double s = (X [j] >= 0) ? 1 : -1 ;
	            if (s != (int) REAL (S [j]))
	            {
	                S [j] = s ;
	                unchanged = FALSE ;
	            }
	        }

	        if (i > 0 && (ainv_norm <= est_old || unchanged == 1))
	        {
	            break ;
	        }

	        for (j = 0 ; j < n ; j++)
	        {
	            X [j] = S [j] ;
	        }

	        /* do a transpose solve */
	        Dklu_tsolve.klu_tsolve (Symbolic, Numeric, n, 1, X, Common) ;

	        /* jnew = the position of the largest entry in X */
	        jnew = 0 ;
	        Xmax = 0 ;
	        for (j = 0 ; j < n ; j++)
	        {
	            /* xj = ABS (X [j]) ;*/
	            ABS (xj, X [j]) ;
	            if (xj > Xmax)
	            {
	                Xmax = xj ;
	                jnew = j ;
	            }
	        }
	        if (i > 0 && jnew == jmax)
	        {
	            /* the position of the largest entry did not change
	             * from the previous iteration */
	            break ;
	        }
	        jmax = jnew ;
	    }

	    /* ---------------------------------------------------------------------- */
	    /* compute another estimate of norm(inv(A),1), and take the largest one */
	    /* ---------------------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
	        CLEAR (X [j]) ;
	        if (j % 2 == 1)
	        {
	            REAL (X [j]) = 1 + ((double) j) / ((double) (n-1)) ;
	        }
	        else
	        {
	            REAL (X [j]) = -1 - ((double) j) / ((double) (n-1)) ;
	        }
	    }

	    Dklu_solve.klu_solve (Symbolic, Numeric, n, 1, (double[]) X, Common) ;

	    est_new = 0.0 ;
	    for (j = 0 ; j < n ; j++)
	    {
	        /* est_new += ABS (X [j]) ;*/
	        ABS (abs_value, X [j]) ;
	        est_new += abs_value ;
	    }
	    est_new = 2 * est_new / (3 * n) ;
	    ainv_norm = MAX (est_new, ainv_norm) ;

	    /* ---------------------------------------------------------------------- */
	    /* compute estimate of condition number */
	    /* ---------------------------------------------------------------------- */

	    Common.condest = ainv_norm * anorm ;
	    return (TRUE) ;
	}

	/**
	 * Compute the flop count for the LU factorization (in Common.flops)
	 *
	 * @param Symbolic
	 * @param Numeric
	 * @param Common
	 * @return TRUE if successful, FALSE otherwise
	 */
	public static int klu_flops(KLU_symbolic Symbolic, KLU_numeric Numeric,
			KLU_common Common)
	{
		double flops = 0 ;
	    int[] R, Ui, Uip, Llen, Ulen ;
	    Unit[][] LUbx ;
	    Unit[] LU ;
	    int k, ulen, p, n, nk, block, nblocks, k1 ;

	    /* ---------------------------------------------------------------------- */
	    /* check inputs */
	    /* ---------------------------------------------------------------------- */

	    if (Common == null)
	    {
	        return (FALSE) ;
	    }
	    Common.flops = EMPTY ;
	    if (Numeric == null || Symbolic == null)
	    {
	        Common.status = KLU_common.KLU_INVALID ;
	        return (FALSE) ;
	    }
	    Common.status = KLU_common.KLU_OK ;

	    /* ---------------------------------------------------------------------- */
	    /* get the contents of the Symbolic object */
	    /* ---------------------------------------------------------------------- */

	    n = Symbolic.n ;
	    R = Symbolic.R ;
	    nblocks = Symbolic.nblocks ;

	    /* ---------------------------------------------------------------------- */
	    /* get the contents of the Numeric object */
	    /* ---------------------------------------------------------------------- */

	    LUbx = (Unit[][]) Numeric.LUbx ;

	    /* ---------------------------------------------------------------------- */
	    /* compute the flop count */
	    /* ---------------------------------------------------------------------- */

	    for (block = 0 ; block < nblocks ; block++)
	    {
	        k1 = R [block] ;
	        nk = R [block+1] - k1 ;
	        if (nk > 1)
	        {
	            Llen = Numeric.Llen + k1 ;
	            Uip  = Numeric.Uip  + k1 ;
	            Ulen = Numeric.Ulen + k1 ;
	            LU = LUbx [block] ;
	            for (k = 0 ; k < nk ; k++)
	            {
	                /* compute kth column of U, and update kth column of A */
	                GET_I_POINTER (LU, Uip, Ui, k) ;
	                ulen = Ulen [k] ;
	                for (p = 0 ; p < ulen ; p++)
	                {
	                    flops += 2 * Llen [Ui [p]] ;
	                }
	                /* gather and divide by pivot to get kth column of L */
	                flops += Llen [k] ;
	            }
	        }
	    }
	    Common.flops = flops ;
	    return (TRUE) ;
	}

	/**
	 * Compute a really cheap estimate of the reciprocal of the condition number,
	 * condition number, min(abs(diag(U))) / max(abs(diag(U))).  If U has a zero
	 * pivot, or a NaN pivot, rcond will be zero.  Takes O(n) time.
	 *
	 * @param Symbolic
	 * @param Numeric
	 * @param Common result in Common.rcond
	 * @return TRUE if successful, FALSE otherwise
	 */
	public static int klu_rcond(KLU_symbolic Symbolic, KLU_numeric Numeric,
		    KLU_common Common)
	{
		double ukk, umin = 0, umax = 0 ;
	    Entry[] Udiag ;
	    int j, n ;

	    /* ---------------------------------------------------------------------- */
	    /* check inputs */
	    /* ---------------------------------------------------------------------- */

	    if (Common == null)
	    {
	        return (FALSE) ;
	    }
	    if (Symbolic == null)
	    {
	        Common.status = KLU_common.KLU_INVALID ;
	        return (FALSE) ;
	    }
	    if (Numeric == null)
	    {
	        Common.rcond = 0 ;
	        Common.status = KLU_common.KLU_SINGULAR ;
	        return (TRUE) ;
	    }
	    Common.status = KLU_common.KLU_OK ;

	    /* ---------------------------------------------------------------------- */
	    /* compute rcond */
	    /* ---------------------------------------------------------------------- */

	    n = Symbolic.n ;
	    Udiag = Numeric.Udiag ;
	    for (j = 0 ; j < n ; j++)
	    {
	        /* get the magnitude of the pivot */
	        ABS (ukk, Udiag [j]) ;
	        if (SCALAR_IS_NAN (ukk) || SCALAR_IS_ZERO (ukk))
	        {
	            /* if NaN, or zero, the rcond is zero */
	            Common.rcond = 0 ;
	            Common.status = KLU_common.KLU_SINGULAR ;
	            return (TRUE) ;
	        }
	        if (j == 0)
	        {
	            /* first pivot entry */
	            umin = ukk ;
	            umax = ukk ;
	        }
	        else
	        {
	            /* subsequent pivots */
	            umin = MIN (umin, ukk) ;
	            umax = MAX (umax, ukk) ;
	        }
	    }

	    Common.rcond = umin / umax ;
	    if (SCALAR_IS_NAN (Common.rcond) || SCALAR_IS_ZERO (Common.rcond))
	    {
	        /* this can occur if umin or umax are Inf or NaN */
	        Common.rcond = 0 ;
	        Common.status = KLU_common.KLU_SINGULAR ;
	    }
	    return (TRUE) ;
	}

}
