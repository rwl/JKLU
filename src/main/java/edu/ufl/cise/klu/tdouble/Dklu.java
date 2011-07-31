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

/**
 * KLU: factorizes P*A into L*U, using the Gilbert-Peierls algorithm[1], with
 * optional symmetric pruning by Eisenstat and Liu[2].  The code is by Tim
 * Davis.  This algorithm is what appears as the default sparse LU routine in
 * MATLAB version 6.0, and still appears in MATLAB 6.5 as[L,U,P] = lu (A).
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
 * Ai[Ap[j] ... Ap[j+1]-1] and the same range of indices in Ax holds the
 * numerical values.  No duplicate entries are allowed.
 *
 * Copyright 2004-2009, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers[1,2].
 *
 * If an input permutation Q is provided, the factorization L*U = A (P,Q)
 * is computed, where P is determined by partial pivoting, and Q is the input
 * ordering.  If the pivot tolerance is less than 1, the "diagonal" entry that
 * KLU attempts to choose is the diagonal of A (Q,Q).  In other words, the
 * input permutation is applied symmetrically to the input matrix.  The output
 * permutation P includes both the partial pivoting ordering and the input
 * permutation.  If Q is null, then it is assumed to be the identity
 * permutation.  Q is not modified.
 *
 *[1] Gilbert, J. R. and Peierls, T., "Sparse Partial Pivoting in Time
 *      Proportional to Arithmetic Operations," SIAM J. Sci. Stat. Comp.,
 *      vol 9, pp.  862-874, 1988.
 *[2] Eisenstat, S. C. and Liu, J. W. H., "Exploiting Structural Symmetry in
 *      Unsymmetric Sparse Symbolic Factorization," SIAM J. Matrix Analysis &
 *      Applic., vol 13, pp.  202-211, 1992.
 */
public class Dklu extends Dklu_internal {

	/**
	 *
	 * @param n A is n-by-n. n must be > 0.
	 * @param Ap size n+1, column pointers for A
	 * @param Ai size nz = Ap[n], row indices for A
	 * @param Ax size nz, values of A
	 * @param Q size n, optional column permutation
	 * @param Lsize estimate of number of nonzeros in L
	 * @param p_LU row indices and values of L and U
	 * @param Udiag size n, diagonal of U
	 * @param Llen size n, column length of L
	 * @param Ulen size n, column length of U
	 * @param Lip size n, column pointers for L
	 * @param Uip size n, column pointers for U
	 * @param P row permutation, size n
	 * @param lnz size of L
	 * @param unz size of U
	 * @param X size n double's, zero on output
	 * @param Work size 5n int's
	 * @param k1 the block of A is from k1 to k2-1
	 * @param PSinv inverse of P from symbolic factorization
	 * @param Rs scale factors for A
	 * @param Offp off-diagonal matrix (modified by this routine)
	 * @param Offi
	 * @param Offx
	 * @param Common
	 * @return
	 */
	public static size_t klu_kernel_factor(int n, int[] Ap, int[] Ai,
			Entry[] Ax, int[] Q, double Lsize,
			Unit[] p_LU, Entry[] Udiag, int[] Llen, int[] Ulen, int[] Lip,
			int[] Uip, int P[], int[] lnz, int[] unz,
			Entry[] X, int[] Work, int k1, int[] PSinv, double[] Rs,
			int[] Offp, int[] Offi, Entry[] Offx, KLU_common Common)
	{
		double maxlnz, dunits;
		Unit[] LU;
		int[] Pinv, Lpend, Stack, Flag, Ap_pos, W;
		int lsize, usize, anz, ok;
		size_t lusize;
		ASSERT (Common != null);

		/* ------------------------------------------------------------------ */
		/* get control parameters, or use defaults */
		/* ------------------------------------------------------------------ */

		n = max(1, n);
		anz = Ap[n+k1] - Ap[k1];

		if (Lsize <= 0)
		{
			Lsize = -Lsize;
			Lsize = max(Lsize, 1.0);
			lsize = Lsize * anz + n;
		}
		else
		{
			lsize = Lsize;
		}

		usize = lsize;

		lsize  = max(n+1, lsize);
		usize  = max(n+1, usize);

		maxlnz = (((double) n) * ((double) n) + ((double) n)) / 2.;
		maxlnz = min(maxlnz, ((double) INT_MAX));
		lsize  = min(maxlnz, lsize);
		usize  = min(maxlnz, usize);

		printf("Welcome to klu: n %d anz %d k1 %d lsize %d usize %d maxlnz %g\n",
			n, anz, k1, lsize, usize, maxlnz);

		/* ------------------------------------------------------------------ */
		/* allocate workspace and outputs */
		/* ------------------------------------------------------------------ */

		/* return arguments are not yet assigned */
		p_LU = null;

		/* these computations are safe from size_t overflow */
		W = Work;
		Pinv = (int[]) W;      W += n;
		Stack = (int[]) W;     W += n;
		Flag = (int[]) W;      W += n;
		Lpend = (int[]) W;     W += n;
		Ap_pos = (int[]) W;    W += n;

		dunits= dunits(Int, lsize) + dunits(Entry, lsize) +
				 dunits(Int, usize) + dunits(Entry, usize);
		lusize = (size_t) dunits;
		ok = !int_overflow(dunits);
		LU = ok != 0 ? Dklu_malloc.klu_malloc(lusize, sizeof(Unit), Common) : null;
		if (LU == null)
		{
			/* out of memory, or problem too large */
			Common.status = KLU_OUT_OF_MEMORY;
			lusize = 0;
			return (lusize);
		}

		/* ------------------------------------------------------------------ */
		/* factorize */
		/* ------------------------------------------------------------------ */

		/* with pruning, and non-recursive depth-first-search */
		lusize = Dklu_kernel.klu_kernel(n, Ap, Ai, Ax, Q, lusize,
				Pinv, P, LU, Udiag, Llen, Ulen, Lip, Uip, lnz, unz,
				X, Stack, Flag, Ap_pos, Lpend,
				k1, PSinv, Rs, Offp, Offi, Offx, Common);

		/* ------------------------------------------------------------------ */
		/* return LU factors, or return nothing if an error occurred */
		/* ------------------------------------------------------------------ */

		if (Common.status < KLU_OK)
		{
			LU = Dklu_free.klu_free(LU, lusize, sizeof(Unit), Common);
			lusize = 0;
		}
		p_LU = LU;
		printf(" in klu noffdiag %d\n", Common.noffdiag);
		return lusize;

	}

	/**
	 * Solve Lx=b.  Assumes L is unit lower triangular and where the unit diagonal
	 * entry is NOT stored.  Overwrites B  with the solution X.  B is n-by-nrhs
	 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
	 * range 1 to 4.
	 *
	 * @param n
	 * @param Lip
	 * @param Llen
	 * @param LU
	 * @param nrhs
	 * @param X right-hand-side on input, solution to Lx=b on output
	 */
	public static void klu_lsolve(int n, int[] Lip, int Llen, Unit[] LU,
			int nrhs, Entry[] X)
	{
		Entry[] x = new Entry[4];
		Entry lik;
		int[] Li;
		Entry[] Lx;
		int k, p, len, i;

		switch (nrhs)
		{

			case 1:
				for (k = 0; k < n; k++)
				{
					x[0] = X[k];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					/* unit diagonal of L is not stored*/
					for (p = 0; p < len; p++)
					{
						/* X[Li[p]] -= Lx[p] * x[0]; */
						mult_sub(X[Li[p]], Lx[p], x[0]);
					}
				}
				break;

			case 2:

				for (k = 0; k < n; k++)
				{
					x[0] = X[2*k    ];
					x[1] = X[2*k + 1];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(X[2*i], lik, x[0]);
						mult_sub(X[2*i + 1], lik, x[1]);
					}
				}
				break;

			case 3:

				for (k = 0; k < n; k++)
				{
					x[0] = X[3*k    ];
					x[1] = X[3*k + 1];
					x[2] = X[3*k + 2];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(X[3*i], lik, x[0]);
						mult_sub(X[3*i + 1], lik, x[1]);
						mult_sub(X[3*i + 2], lik, x[2]);
					}
				}
				break;

			case 4:

				for (k = 0; k < n; k++)
				{
					x[0] = X[4*k    ];
					x[1] = X[4*k + 1];
					x[2] = X[4*k + 2];
					x[3] = X[4*k + 3];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(X[4*i], lik, x[0]);
						mult_sub(X[4*i + 1], lik, x[1]);
						mult_sub(X[4*i + 2], lik, x[2]);
						mult_sub(X[4*i + 3], lik, x[3]);
					}
				}
				break;

		}

	}

	/**
	 * Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
	 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
	 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
	 * range 1 to 4.
	 *
	 * @param n
	 * @param Uip
	 * @param Ulen
	 * @param LU
	 * @param Udiag
	 * @param nrhs
	 * @param X right-hand-side on input, solution to Ux=b on output
	 */
	public static void klu_usolve(int n, int[] Uip, int[] Ulen, Unit[] LU,
			Entry[] Udiag, int nrhs, Entry[] X)
	{
		Entry[] x = new Entry[4];
		Entry uik, ukk;
		int[] Ui;
		Entry[] Ux;
		int k, p, len, i;

		switch (nrhs)
		{

			case 1:

				for (k = n-1; k >= 0; k--)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					/* x[0] = X[k] / Udiag[k]; */
					div(x[0], X[k], Udiag[k]);
					X[k] = x[0];
					for (p = 0; p < len; p++)
					{
						/* X[Ui[p]] -= Ux[p] * x[0]; */
						mult_sub(X[Ui[p]], Ux[p], x[0]);

					}
				}

				break;

			case 2:

				for (k = n-1; k >= 0; k--)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					ukk = Udiag[k];
					/* x[0] = X[2*k    ] / ukk;
					x[1] = X[2*k + 1] / ukk; */
					div(x[0], X[2*k], ukk);
					div(x[1], X[2*k + 1], ukk);

					X[2*k    ] = x[0];
					X[2*k + 1] = x[1];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];
						/* X[2*i    ] -= uik * x[0];
						X[2*i + 1] -= uik * x[1]; */
						mult_sub(X[2*i], uik, x[0]);
						mult_sub(X[2*i + 1], uik, x[1]);
					}
				}

				break;

			case 3:

				for (k = n-1; k >= 0; k--)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					ukk = Udiag[k];

					div(x[0], X[3*k], ukk);
					div(x[1], X[3*k + 1], ukk);
					div(x[2], X[3*k + 2], ukk);

					X[3*k    ] = x[0];
					X[3*k + 1] = x[1];
					X[3*k + 2] = x[2];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];
						mult_sub(X[3*i], uik, x[0]);
						mult_sub(X[3*i + 1], uik, x[1]);
						mult_sub(X[3*i + 2], uik, x[2]);
					}
				}

				break;

			case 4:

				for (k = n-1; k >= 0; k--)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					ukk = Udiag[k];

					div(x[0], X[4*k], ukk);
					div(x[1], X[4*k + 1], ukk);
					div(x[2], X[4*k + 2], ukk);
					div(x[3], X[4*k + 3], ukk);

					X[4*k    ] = x[0];
					X[4*k + 1] = x[1];
					X[4*k + 2] = x[2];
					X[4*k + 3] = x[3];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];

						mult_sub(X[4*i], uik, x[0]);
						mult_sub(X[4*i + 1], uik, x[1]);
						mult_sub(X[4*i + 2], uik, x[2]);
						mult_sub(X[4*i + 3], uik, x[3]);
					}
				}

				break;

		}

	}

	/**
	 * Solve L'x=b.  Assumes L is unit lower triangular and where the unit diagonal
	 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
	 * and is stored in ROW form with row dimension nrhs.  nrhs must in the
	 * range 1 to 4.
	 *
	 * @param n
	 * @param Lip
	 * @param Llen
	 * @param LU
	 * @param nrhs
	 * @param X right-hand-side on input, solution to L'x=b on output
	 */
	public static void klu_ltsolve(int n, int[] Lip, int[] Llen, Unit[] LU,
			int nrhs, Entry[] X)
	{
		Entry[] x = Entry[4];
		Entry lik;
		int[] Li;
		Entry[] Lx;
		int k, p, len, i;

		switch (nrhs)
		{

			case 1:

				for (k = n-1; k >= 0; k--)
				{
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					x[0] = X[k];
					for (p = 0; p < len; p++)
					{
						/*x[0] -= Lx[p] * X[Li[p]];*/
						mult_sub(x[0], Lx[p], X[Li[p]]);
					}
					X[k] = x[0];
				}
				break;

			case 2:

				for (k = n-1; k >= 0; k--)
				{
					x[0] = X[2*k    ];
					x[1] = X[2*k + 1];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(x[0], lik, X[2*i]);
						mult_sub(x[1], lik, X[2*i + 1]);
					}
					X[2*k    ] = x[0];
					X[2*k + 1] = x[1];
				}
				break;

			case 3:

				for (k = n-1; k >= 0; k--)
				{
					x[0] = X[3*k    ];
					x[1] = X[3*k + 1];
					x[2] = X[3*k + 2];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(x[0], lik, X[3*i]);
						mult_sub(x[1], lik, X[3*i + 1]);
						mult_sub(x[2], lik, X[3*i + 2]);
					}
					X[3*k    ] = x[0];
					X[3*k + 1] = x[1];
					X[3*k + 2] = x[2];
				}
				break;

			case 4:

				for (k = n-1; k >= 0; k--)
				{
					x[0] = X[4*k    ];
					x[1] = X[4*k + 1];
					x[2] = X[4*k + 2];
					x[3] = X[4*k + 3];
					get_pointer(LU, Lip, Llen, Li, Lx, k, len);
					for (p = 0; p < len; p++)
					{
						i = Li[p];
						lik = Lx[p];
						mult_sub(x[0], lik, X[4*i]);
						mult_sub(x[1], lik, X[4*i + 1]);
						mult_sub(x[2], lik, X[4*i + 2]);
						mult_sub(x[3], lik, X[4*i + 3]);
					}
					X[4*k    ] = x[0];
					X[4*k + 1] = x[1];
					X[4*k + 2] = x[2];
					X[4*k + 3] = x[3];
				}
				break;
		}
	}

	/**
	 * Solve U'x=b.  Assumes U is non-unit upper triangular and where the diagonal
	 * entry is stored (and appears last in each column of U).  Overwrites B
	 * with the solution X.  B is n-by-nrhs and is stored in ROW form with row
	 * dimension nrhs.  nrhs must be in the range 1 to 4.
	 *
	 * @param n
	 * @param Uip
	 * @param Ulen
	 * @param LU
	 * @param Udiag
	 * @param nrhs
	 * @param X right-hand-side on input, solution to Ux=b on output
	 */
	public static void klu_utsolve(int n, int[] Uip, int[] Ulen, Unit[] LU,
			Entry[] Udiag, int nrhs, Entry[] X)
	{
		Entry[] x = new Entry[4];
		Entry uik, ukk;
		int k, p, len, i;
		int[] Ui;
		Entry[] Ux;

		switch (nrhs)
		{

			case 1:

				for (k = 0; k < n; k++)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					x[0] = X[k];
					for (p = 0; p < len; p++)
					{
						/* x[0] -= Ux[p] * X[Ui[p]]; */
						mult_sub(x[0], Ux[p], X[Ui[p]]);
					}
					ukk = Udiag[k];
					div(X[k], x[0], ukk);
				}
				break;

			case 2:

				for (k = 0; k < n; k++)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					x[0] = X[2*k    ];
					x[1] = X[2*k + 1];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];
						mult_sub(x[0], uik, X[2*i]);
						mult_sub(x[1], uik, X[2*i + 1]);
					}
					ukk = Udiag[k];
					div(X[2*k], x[0], ukk);
					div(X[2*k + 1], x[1], ukk);
				}
				break;

			case 3:

				for (k = 0; k < n; k++)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					x[0] = X[3*k    ];
					x[1] = X[3*k + 1];
					x[2] = X[3*k + 2];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];
						mult_sub(x[0], uik, X[3*i]);
						mult_sub(x[1], uik, X[3*i + 1]);
						mult_sub(x[2], uik, X[3*i + 2]);
					}
					ukk = Udiag[k];
					div(X[3*k], x[0], ukk);
					div(X[3*k + 1], x[1], ukk);
					div(X[3*k + 2], x[2], ukk);
				}
				break;

			case 4:

				for (k = 0; k < n; k++)
				{
					get_pointer(LU, Uip, Ulen, Ui, Ux, k, len);
					x[0] = X[4*k    ];
					x[1] = X[4*k + 1];
					x[2] = X[4*k + 2];
					x[3] = X[4*k + 3];
					for (p = 0; p < len; p++)
					{
						i = Ui[p];
						uik = Ux[p];
						mult_sub(x[0], uik, X[4*i]);
						mult_sub(x[1], uik, X[4*i + 1]);
						mult_sub(x[2], uik, X[4*i + 2]);
						mult_sub(x[3], uik, X[4*i + 3]);
					}
					ukk = Udiag[k];
					div(X[4*k], x[0], ukk);
					div(X[4*k + 1], x[1], ukk);
					div(X[4*k + 2], x[2], ukk);
					div(X[4*k + 3], x[3], ukk);
				}
				break;
		}
	}

}
