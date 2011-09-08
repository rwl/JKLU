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

package edu.ufl.cise.klu.common;

/**
 * Numeric object - contains the factors computed by klu_factor.
 *
 * 64-bit version.
 */
public class KLU_l_numeric
{

	/* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    long n;             /* A is n-by-n */
    long nblocks;       /* number of diagonal blocks */
    long lnz;           /* actual nz in L, including diagonal */
    long unz;           /* actual nz in U, including diagonal */
    long max_lnz_block; /* max actual nz in L in any one block, incl. diag */
    long max_unz_block; /* max actual nz in U in any one block, incl. diag */
    long[] Pnum;        /* size n. final pivot permutation */
    long[] Pinv;        /* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    long[] Lip;         /* size n. pointers into LUbx[block] for L */
    long[] Uip;         /* size n. pointgers into LUbx[block] for U */
    int[] Llen;        /* size n. Llen [k] = # of entries in kth column of L */
    long[] Ulen;        /* size n. Ulen [k] = # of entries in kth column of U */
    double[][] LUbx;      /* L and U indices and entries (excl. diagonal of U) */
    long[] LUsize;    /* size of each LUbx [block], in sizeof (Unit) */
    double[] Udiag;       /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double[] Rs;        /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    long[] worksize;  /* size (in bytes) of Work */
    double[] Work;        /* workspace */
    double[] Xwork;       /* alias into Numeric->Work */
    long[] Iwork;       /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    long[] Offp;        /* size n+1, column pointers */
    long[] Offi;        /* size nzoff, row indices */
    double[] Offx;        /* size nzoff, numerical values */
    long nzoff;

}
