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
 * Numeric object - contains the factors computed by klu_factor.
 */
public class KLU_numeric
{

	/* LU factors of each block, the pivot row permutation, and the
	 * entries in the off-diagonal blocks */

	public int n;             /* A is n-by-n */
	public int nblocks;       /* number of diagonal blocks */
	public int lnz;           /* actual nz in L, including diagonal */
	public int unz;           /* actual nz in U, including diagonal */
	public int max_lnz_block; /* max actual nz in L in any one block, incl. diag */
	public int max_unz_block; /* max actual nz in U in any one block, incl. diag */
	public int[] Pnum;        /* size n. final pivot permutation */
	public int[] Pinv;        /* size n. inverse of final pivot permutation */

	/* LU factors of each block */
	public int[] Lip;         /* size n. pointers into LUbx[block] for L */
	public int[] Uip;         /* size n. pointers into LUbx[block] for U */
	public int[] Llen;        /* size n. Llen [k] = # of entries in kth column of L */
	public int[] Ulen;        /* size n. Ulen [k] = # of entries in kth column of U */
	public double[] LUbx;     /* L and U indices and entries (excl. diagonal of U) */
	public int[] LUsize;   /* size of each LUbx [block], in sizeof (Unit) */
	public double[] Udiag;      /* diagonal of U */

	/* scale factors; can be NULL if no scaling */
	public double[] Rs;       /* size n. Rs [i] is scale factor for row i */

	/* permanent workspace for factorization and solve */
	public int worksize; /* size (in bytes) of Work */
	public double[] Work;       /* workspace */
	public double[] Xwork;      /* alias into Numeric->Work */
	public int[] Iwork;       /* alias into Numeric->Work */

	/* off-diagonal entries in a conventional compressed-column sparse matrix */
	public int[] Offp;        /* size n+1, column pointers */
	public int[] Offi;        /* size nzoff, row indices */
	public double[] Offx;       /* size nzoff, numerical values */
	public int nzoff;

}
