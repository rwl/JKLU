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

import static edu.ufl.cise.klu.tdouble.Dklu_dump.klu_valid;

/**
 * Solve Ax=b using the symbolic and numeric objects from KLU_analyze
 * (or KLU_analyze_given) and KLU_factor.  Note that no iterative refinement is
 * performed.  Uses Numeric.Xwork as workspace (undefined on input and output),
 * of size 4n double's (note that columns 2 to 4 of Xwork overlap with
 * Numeric.Iwork).
 */
public class Dklu_solve extends Dklu_internal {

	/**
	 *
	 * @param Symbolic
	 * @param Numeric
	 * @param d leading dimension of B
	 * @param nrhs number of right-hand-sides
	 * @param B right-hand-side on input, overwritten with solution to Ax=b on
	 * output. Size n*nrhs, in column-oriented form, with leading dimension d.
	 * @param Common
	 * @return
	 */
	public static int klu_solve(KLU_symbolic Symbolic, KLU_numeric Numeric,
			int d, int nrhs, double[] B, KLU_common Common)
	{
		double offik, s ;
		double[] x = new double[4] ;
		double rs ;
		double[] Offx, X, Bz, Udiag, Rs ;
		int[] Q, R, Pnum, Offp, Offi, Lip, Uip, Llen, Ulen ;
		double[][] LUbx ;
		int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i ;

		/* ---------------------------------------------------------------------- */
		/* check inputs */
		/* ---------------------------------------------------------------------- */

		if (Common == null)
		{
			return (FALSE) ;
		}
		if (Numeric == null || Symbolic == null || d < Symbolic.n || nrhs < 0 ||
			B == null)
		{
			Common.status = KLU_INVALID ;
			return (FALSE) ;
		}
		Common.status = KLU_OK ;

		/* ---------------------------------------------------------------------- */
		/* get the contents of the Symbolic object */
		/* ---------------------------------------------------------------------- */

		Bz = B ;
		n = Symbolic.n ;
		nblocks = Symbolic.nblocks ;
		Q = Symbolic.Q ;
		R = Symbolic.R ;

		/* ---------------------------------------------------------------------- */
		/* get the contents of the Numeric object */
		/* ---------------------------------------------------------------------- */

		ASSERT (nblocks == Numeric.nblocks) ;
		Pnum = Numeric.Pnum ;
		Offp = Numeric.Offp ;
		Offi = Numeric.Offi ;
		Offx = Numeric.Offx ;

		Lip  = Numeric.Lip ;
		Llen = Numeric.Llen ;
		Uip  = Numeric.Uip ;
		Ulen = Numeric.Ulen ;
		LUbx = Numeric.LUbx ;
		Udiag = Numeric.Udiag ;

		Rs = Numeric.Rs ;
		X = Numeric.Xwork ;

		ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

		/* ---------------------------------------------------------------------- */
		/* solve in chunks of 4 columns at a time */
		/* ---------------------------------------------------------------------- */

		for (chunk = 0 ; chunk < nrhs ; chunk += 4)
		{

			/* ------------------------------------------------------------------ */
			/* get the size of the current chunk */
			/* ------------------------------------------------------------------ */

			nr = MIN (nrhs - chunk, 4) ;

			/* ------------------------------------------------------------------ */
			/* scale and permute the right hand side, X = P*(R\B) */
			/* ------------------------------------------------------------------ */

			if (Rs == null)
			{

				/* no scaling */
				switch (nr)
				{

					case 1:

						for (k = 0 ; k < n ; k++)
						{
							X [k] = Bz [Pnum [k]] ;
						}
						break ;

					case 2:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							X [2*k    ] = Bz [i      ] ;
							X [2*k + 1] = Bz  [i + d  ] ;
						}
						break ;

					case 3:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							X [3*k    ] = Bz [i      ] ;
							X [3*k + 1] = Bz [i + d  ] ;
							X [3*k + 2] = Bz [i + d*2] ;
						}
						break ;

					case 4:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							X [4*k    ] = Bz [i      ] ;
							X [4*k + 1] = Bz [i + d  ] ;
							X [4*k + 2] = Bz [i + d*2] ;
							X [4*k + 3] = Bz [i + d*3] ;
						}
						break ;
				}

			}
			else
			{

				switch (nr)
				{

					case 1:

						for (k = 0 ; k < n ; k++)
						{
							X [k] = Bz  [Pnum [k]] / Rs [k] ;
							//SCALE_DIV_ASSIGN (X [k], Bz  [Pnum [k]], Rs [k]) ;
						}
						break ;

					case 2:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							rs = Rs [k] ;
							X [2*k] =  Bz [i] / rs ;
							//SCALE_DIV_ASSIGN (X [2*k], Bz [i], rs) ;
							X [2*k + 1] = Bz [i + d] / rs ;
							//SCALE_DIV_ASSIGN (X [2*k + 1], Bz [i + d], rs) ;
						}
						break ;

					case 3:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							rs = Rs [k] ;
							X [3*k] = Bz [i] / rs ;
							//SCALE_DIV_ASSIGN (X [3*k], Bz [i], rs) ;
							X [3*k + 1] = Bz [i + d] / rs ;
							//SCALE_DIV_ASSIGN (X [3*k + 1], Bz [i + d], rs) ;
							X [3*k + 2] = Bz [i + d*2] / rs ;
							//SCALE_DIV_ASSIGN (X [3*k + 2], Bz [i + d*2], rs) ;
						}
						break ;

					case 4:

						for (k = 0 ; k < n ; k++)
						{
							i = Pnum [k] ;
							rs = Rs [k] ;
							X [4*k] = Bz [i] / rs ;
							//SCALE_DIV_ASSIGN (X [4*k], Bz [i], rs) ;
							X [4*k + 1] = Bz [i + d] / rs ;
							//SCALE_DIV_ASSIGN (X [4*k + 1], Bz [i + d], rs) ;
							X [4*k + 2] = Bz [i + d*2] / rs ;
							//SCALE_DIV_ASSIGN (X [4*k + 2], Bz [i + d*2], rs) ;
							X [4*k + 3] = Bz [i + d*3] / rs ;
							//SCALE_DIV_ASSIGN (X [4*k + 3], Bz [i + d*3], rs) ;
						}
						break ;
				}
			}

			/* ------------------------------------------------------------------ */
			/* solve X = (L*U + Off)\X */
			/* ------------------------------------------------------------------ */

			for (block = nblocks-1 ; block >= 0 ; block--)
			{

				/* -------------------------------------------------------------- */
				/* the block of size nk is from rows/columns k1 to k2-1 */
				/* -------------------------------------------------------------- */

				k1 = R [block] ;
				k2 = R [block+1] ;
				nk = k2 - k1 ;
				PRINTF ("solve %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk) ;

				/* solve the block system */
				if (nk == 1)
				{
					s = Udiag [k1] ;
					switch (nr)
					{

						case 1:
							X [k1] = X [k1] / s ;
							//DIV (X [k1], X [k1], s) ;
							break ;

						case 2:
							X [2*k1] = X [2*k1] / s ;
							//DIV (X [2*k1], X [2*k1], s) ;
							X [2*k1 + 1] = X [2*k1 + 1] / s ;
							//DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
							break ;

						case 3:
							X [3*k1] = X [3*k1] / s ;
							//DIV (X [3*k1], X [3*k1], s) ;
							X [3*k1 + 1] = X [3*k1 + 1] / s ;
							//DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
							X [3*k1 + 2] = X [3*k1 + 2] / s ;
							//DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
							break ;

						case 4:
							X [4*k1] = X [4*k1] / s ;
							//DIV (X [4*k1], X [4*k1], s) ;
							X [4*k1 + 1] = X [4*k1 + 1] / s ;
							//DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
							X [4*k1 + 2] = X [4*k1 + 2] / s ;
							//DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
							X [4*k1 + 3] = X [4*k1 + 3] / s ;
							//DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
							break ;

					}
				}
				else
				{
					Dklu.klu_lsolve (nk, Lip + k1, Llen + k1,
							LUbx [block], nr, X + nr*k1) ;
					Dklu.klu_usolve (nk, Uip + k1, Ulen + k1,
							LUbx [block], Udiag + k1, nr, X + nr*k1) ;
				}

				/* -------------------------------------------------------------- */
				/* block back-substitution for the off-diagonal-block entries */
				/* -------------------------------------------------------------- */

				if (block > 0)
				{
					switch (nr)
					{

						case 1:

							for (k = k1 ; k < k2 ; k++)
							{
								pend = Offp [k+1] ;
								x [0] = X [k] ;
								for (p = Offp [k] ; p < pend ; p++)
								{
									X [Offi [p]] -= Offx [p] * x [0] ;
									//MULT_SUB (X [Offi [p]], Offx [p], x [0]) ;
								}
							}
							break ;

						case 2:

							for (k = k1 ; k < k2 ; k++)
							{
								pend = Offp [k+1] ;
								x [0] = X [2*k    ] ;
								x [1] = X [2*k + 1] ;
								for (p = Offp [k] ; p < pend ; p++)
								{
									i = Offi [p] ;
									offik = Offx [p] ;
									X [2*i] -= offik * x [0] ;
									//MULT_SUB (X [2*i], offik, x [0]) ;
									X [2*i + 1] -= offik * x [1] ;
									//MULT_SUB (X [2*i + 1], offik, x [1]) ;
								}
							}
							break ;

						case 3:

							for (k = k1 ; k < k2 ; k++)
							{
								pend = Offp [k+1] ;
								x [0] = X [3*k    ] ;
								x [1] = X [3*k + 1] ;
								x [2] = X [3*k + 2] ;
								for (p = Offp [k] ; p < pend ; p++)
								{
									i = Offi [p] ;
									offik = Offx [p] ;
									X [3*i] -= offik * x [0] ;
									//MULT_SUB (X [3*i], offik, x [0]) ;
									X [3*i + 1] -= offik * x [1] ;
									//MULT_SUB (X [3*i + 1], offik, x [1]) ;
									X [3*i + 2] -= offik * x [2] ;
									//MULT_SUB (X [3*i + 2], offik, x [2]) ;
								}
							}
							break ;

						case 4:

							for (k = k1 ; k < k2 ; k++)
							{
								pend = Offp [k+1] ;
								x [0] = X [4*k    ] ;
								x [1] = X [4*k + 1] ;
								x [2] = X [4*k + 2] ;
								x [3] = X [4*k + 3] ;
								for (p = Offp [k] ; p < pend ; p++)
								{
									i = Offi [p] ;
									offik = Offx [p] ;
									X [4*i] -= offik * x [0] ;
									//MULT_SUB (X [4*i], offik, x [0]) ;
									X [4*i + 1] -= offik * x [1] ;
									//MULT_SUB (X [4*i + 1], offik, x [1]) ;
									X [4*i + 2] -= offik * x [2] ;
									//MULT_SUB (X [4*i + 2], offik, x [2]) ;
									X [4*i + 3] -= offik * x [3] ;
									//MULT_SUB (X [4*i + 3], offik, x [3]) ;
								}
							}
							break ;
					}
				}
			}

			/* ------------------------------------------------------------------ */
			/* permute the result, Bz  = Q*X */
			/* ------------------------------------------------------------------ */

			switch (nr)
			{

				case 1:

					for (k = 0 ; k < n ; k++)
					{
						Bz  [Q [k]] = X [k] ;
					}
					break ;

				case 2:

					for (k = 0 ; k < n ; k++)
					{
						i = Q [k] ;
						Bz  [i      ] = X [2*k    ] ;
						Bz  [i + d  ] = X [2*k + 1] ;
					}
					break ;

				case 3:

					for (k = 0 ; k < n ; k++)
					{
						i = Q [k] ;
						Bz  [i      ] = X [3*k    ] ;
						Bz  [i + d  ] = X [3*k + 1] ;
						Bz  [i + d*2] = X [3*k + 2] ;
					}
					break ;

				case 4:

					for (k = 0 ; k < n ; k++)
					{
						i = Q [k] ;
						Bz  [i      ] = X [4*k    ] ;
						Bz  [i + d  ] = X [4*k + 1] ;
						Bz  [i + d*2] = X [4*k + 2] ;
						Bz  [i + d*3] = X [4*k + 3] ;
					}
					break ;
			}

			/* ------------------------------------------------------------------ */
			/* go to the next chunk of B */
			/* ------------------------------------------------------------------ */

			Bz += d*4 ;
		}
		return (TRUE) ;
	}

}
