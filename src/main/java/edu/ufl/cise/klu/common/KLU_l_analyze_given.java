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
 * Analyzes a matrix using given P and Q.
 *
 * 64-bit version.
 */
public class KLU_l_analyze_given
{

	/* Order the matrix with BTF (or not), then use natural or given ordering
	 * P and Q on the blocks.  P and Q are interpretted as identity
	 * if NULL. */

	/**
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param P size n, user's row permutation (may be NULL)
	 * @param Q size n, user's column permutation (may be NULL)
	 * @param Common
	 */
	public static KLU_l_symbolic klu_l_analyze_given(long n, long[] Ap,
			long[] Ai, long[] P, long[] Q, KLU_l_common Common)
	{
		return null;
	}

}
