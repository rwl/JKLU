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
 * KLU control parameters and statistics.
 *
 * 64-bit version.
 */
public class KLU_l_common
{

	double tol, memgrow, initmem_amd, initmem, maxwork;
    long btf, ordering, scale;
//    void *(*malloc_memory) (size_t);
//    void *(*realloc_memory) (void *, size_t);
//    void (*free_memory) (void *);
//    void *(*calloc_memory) (size_t, size_t);
//    UF_long (*user_order) (UF_long, UF_long *, UF_long *, UF_long *,
//        struct klu_l_common_struct *);
    Object user_data;
    long halt_if_singular;
    long status, nrealloc, structural_rank, numerical_rank, singular_col,
        noffdiag;
    double flops, rcond, condest, rgrowth, work;
    long memusage, mempeak;

}
