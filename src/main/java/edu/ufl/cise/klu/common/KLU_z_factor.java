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
 * Factors a matrix using the klu_analyze results.
 *
 * Complex version.
 */
public class KLU_z_factor
{

	/**
	 *
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz, row indices
	 * @param Ax size 2*nz, numerical values (real,imag pairs)
	 * @param Symbolic
	 * @param Common
	 * @return KLU_OK if OK, < 0 if error
	 */
	public static KLU_numeric klu_factor(int[] Ap, int[] Ai, double[] Ax,
		    KLU_symbolic Symbolic, KLU_common Common)
	{
		return null;
	}

}
