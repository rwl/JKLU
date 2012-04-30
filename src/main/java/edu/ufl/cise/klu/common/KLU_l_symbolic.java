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

package edu.ufl.cise.klu.common;

/**
 * Symbolic object - contains the pre-ordering computed by klu_analyze.
 *
 * 64-bit version.
 */
public class KLU_l_symbolic
{
    /* A (P,Q) is in upper block triangular form.  The kth block goes from
     * row/col index R [k] to R [k+1]-1.  The estimated number of nonzeros
     * in the L factor of the kth block is Lnz [k].
     */

    /* only computed if the AMD ordering is chosen: */
    double symmetry;   /* symmetry of largest block */
    double est_flops;  /* est. factorization flop count */
    double lnz, unz;   /* estimated nz in L and U, including diagonals */
    double[] Lnz;      /* size n, but only Lnz [0..nblocks-1] is used */

    /* computed for all orderings: */
    long
        n,              /* input matrix A is n-by-n */
        nz,             /* # entries in input matrix */
        nzoff,          /* nz in off-diagonal blocks */
        nblocks,        /* number of blocks */
        maxblock,       /* size of largest block */
        ordering,       /* ordering used (AMD, COLAMD, or GIVEN) */
        do_btf;         /* whether or not BTF preordering was requested */

    long[]
        P,              /* size n */
        Q,              /* size n */
        R;              /* size n+1, but only R [0..nblocks] is used */

    /* only computed if BTF preordering requested */
    int structural_rank;   /* 0 to n-1 if the matrix is structurally rank
                        * deficient.  -1 if not computed.  n if the matrix has
                        * full structural rank */

}
