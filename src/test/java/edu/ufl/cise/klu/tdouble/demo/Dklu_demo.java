package edu.ufl.cise.klu.tdouble.demo;

/**
 * Read in a Matrix Market matrix(using CHOLMOD) and solve a linear system.
 */
public class Dklu_demo {

	private static double real(double[] X, int i)
	{
		return X[2*i];
	}

	private static double imag(double[] X, int i)
	{
		return X[2*i + 1];
	}

	private static double cabs(double[] X, int i)
	{
		return Math.sqrt(real(X, i) * real(X, i) + imag(X, i) * imag(X, i));
	}

	private static double max(double a, double b)
	{
		return (a) > (b) ? (a) : (b);
	}

	/**
	 *
	 * @param n A is n-by-n
	 * @param Ap size n+1, column pointers
	 * @param Ai size nz = Ap [n], row indices
	 * @param Ax size nz, numerical values
	 * @param isreal nonzero if A is real, 0 otherwise
	 * @param B size n, right-hand-side
	 * @param X size n, solution to Ax=b
	 * @param R size n, residual r = b-A*x
	 * @param lunz nnz(L+U+F)
	 * @param rnorm norm(b-A*x,1) / norm(A,1)
	 * @param Common default parameters and statistics
	 * @return 1 if successful, 0 otherwise
	 */
	public static int klu_backslash(int n, int[] Ap, int[] Ai, double[] Ax,
			int isreal, double[] B, double[] X, double[] R, int lunz,
			double rnorm, Dklu_common Common)
	{
		double anorm = 0, asum;
		klu_symbolic *Symbolic;
		klu_numeric *Numeric;
		int i, j, p;

		if (!Ap || !Ai || !Ax || !B || !X || !B) return(0);

		/* ---------------------------------------------------------------------- */
		/* symbolic ordering and analysis */
		/* ---------------------------------------------------------------------- */

		Symbolic = klu_analyze(n, Ap, Ai, Common);
		if (!Symbolic) return(0);

		if (isreal)
		{

			/* ------------------------------------------------------------------ */
			/* factorization */
			/* ------------------------------------------------------------------ */

			Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common);
			if (!Numeric)
			{
				klu_free_symbolic(&Symbolic, Common);
				return(0);
			}

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

			klu_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common);
			klu_condest(Ap, Ax, Symbolic, Numeric, Common);
			klu_rcond(Symbolic, Numeric, Common);
			klu_flops(Symbolic, Numeric, Common);
			*lunz = Numeric.lnz + Numeric.unz - n +
				((Numeric.Offp) ?(Numeric.Offp [n]) : 0);

			/* ------------------------------------------------------------------ */
			/* solve Ax=b */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < n; i++)
			{
				X [i] = B [i];
			}
			klu_solve(Symbolic, Numeric, n, 1, X, Common);

			/* ------------------------------------------------------------------ */
			/* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < n; i++)
			{
				R [i] = B [i];
			}
			for(j = 0; j < n; j++)
			{
				asum = 0;
				for(p = Ap [j]; p < Ap [j+1]; p++)
				{
					/* R(i) -= A(i,j) * X(j) */
					R [Ai [p]] -= Ax [p] * X [j];
					asum += fabs(Ax [p]);
				}
				anorm = MAX(anorm, asum);
			}
			*rnorm = 0;
			for(i = 0; i < n; i++)
			{
				*rnorm = MAX(*rnorm, fabs(R [i]));
			}

			/* ------------------------------------------------------------------ */
			/* free numeric factorization */
			/* ------------------------------------------------------------------ */

			klu_free_numeric(&Numeric, Common);

		}
		else
		{

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

			Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
			if (!Numeric)
			{
				klu_free_symbolic(&Symbolic, Common);
				return(0);
			}

			/* ------------------------------------------------------------------ */
			/* statistics */
			/* ------------------------------------------------------------------ */

			klu_z_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common);
			klu_z_condest(Ap, Ax, Symbolic, Numeric, Common);
			klu_z_rcond(Symbolic, Numeric, Common);
			klu_z_flops(Symbolic, Numeric, Common);
			*lunz = Numeric.lnz + Numeric.unz - n +
				((Numeric.Offp) ?(Numeric.Offp [n]) : 0);

			/* ------------------------------------------------------------------ */
			/* solve Ax=b */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < 2*n; i++)
			{
				X [i] = B [i];
			}
			klu_z_solve(Symbolic, Numeric, n, 1, X, Common);

			/* ------------------------------------------------------------------ */
			/* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < 2*n; i++)
			{
				R [i] = B [i];
			}
			for(j = 0; j < n; j++)
			{
				asum = 0;
				for(p = Ap [j]; p < Ap [j+1]; p++)
				{
					/* R(i) -= A(i,j) * X(j) */
					i = Ai [p];
					REAL(R,i) -= REAL(Ax,p) * REAL(X,j) - IMAG(Ax,p) * IMAG(X,j);
					IMAG(R,i) -= IMAG(Ax,p) * REAL(X,j) + REAL(Ax,p) * IMAG(X,j);
					asum += CABS(Ax, p);
				}
				anorm = MAX(anorm, asum);
			}
			*rnorm = 0;
			for(i = 0; i < n; i++)
			{
				*rnorm = MAX(*rnorm, CABS(R, i));
			}

			/* ------------------------------------------------------------------ */
			/* free numeric factorization */
			/* ------------------------------------------------------------------ */

			klu_z_free_numeric(&Numeric, Common);
		}

		/* ---------------------------------------------------------------------- */
		/* free symbolic analysis, and residual */
		/* ---------------------------------------------------------------------- */

		klu_free_symbolic(&Symbolic, Common);
		return(1);
	}

	/**
	 * Given a sparse matrix A, set up a right-hand-side and solve X = A\b.
	 */
	public static void klu_demo(int n, int[] Ap, int[] Ai, double[] Ax,
			int isreal)
	{
		double rnorm;
		klu_common Common;
		double *B, *X, *R;
		int i, lunz;

		System.out.printf("KLU: %s, version: %d.%d.%d\n", KLU_DATE, KLU_MAIN_VERSION,
			KLU_SUB_VERSION, KLU_SUBSUB_VERSION);

		/* ---------------------------------------------------------------------- */
		/* set defaults */
		/* ---------------------------------------------------------------------- */

		klu_defaults(Common);

		/* ---------------------------------------------------------------------- */
		/* create a right-hand-side */
		/* ---------------------------------------------------------------------- */

		if (isreal)
		{
			/* B = 1 +(1:n)/n */
			B = klu_malloc(n, sizeof(Double), Common);
			X = klu_malloc(n, sizeof(Double), Common);
			R = klu_malloc(n, sizeof(Double), Common);
			if (B)
			{
				for(i = 0; i < n; i++)
				{
					B [i] = 1 +((double) i+1) /((double) n);
				}
			}
		}
		else
		{
			/* real(B) = 1 +(1:n)/n, imag(B) = (n:-1:1)/n */
			B = klu_malloc(n, 2 * sizeof(Double), Common);
			X = klu_malloc(n, 2 * sizeof(Double), Common);
			R = klu_malloc(n, 2 * sizeof(Double), Common);
			if (B)
			{
				for(i = 0; i < n; i++)
				{
					REAL(B, i) = 1 +((double) i+1) /((double) n);
					IMAG(B, i) = ((double) n-i) /((double) n);
				}
			}
		}

		/* ---------------------------------------------------------------------- */
		/* X = A\b using KLU and print statistics */
		/* ---------------------------------------------------------------------- */

		if (!klu_backslash(n, Ap, Ai, Ax, isreal, B, X, R, lunz, rnorm, Common))
		{
			System.out.printf("KLU failed\n");
		}
		else
		{
			System.out.printf("n %d nnz(A) %d nnz(L+U+F) %d resid %g\n" +
				"recip growth %g condest %g rcond %g flops %g\n",
				n, Ap [n], lunz, rnorm, Common.rgrowth, Common.condest,
				Common.rcond, Common.flops);
		}

		/* ---------------------------------------------------------------------- */
		/* free the problem */
		/* ---------------------------------------------------------------------- */

		if (isreal)
		{
			klu_free(B, n, sizeof(Double), Common);
			klu_free(X, n, sizeof(Double), Common);
			klu_free(R, n, sizeof(Double), Common);
		}
		else
		{
			klu_free(B, 2*n, sizeof(Double), Common);
			klu_free(X, 2*n, sizeof(Double), Common);
			klu_free(R, 2*n, sizeof(Double), Common);
		}
		System.out.printf("peak memory usage: %g bytes\n\n",(double)(Common.mempeak));
	}

	/**
	 * Read in a sparse matrix in Matrix Market format using CHOLMOD, and then
	 * solve Ax=b with KLU.  Note that CHOLMOD is only used to read the matrix.
	 */
	public static void main(String[] args) {
		cholmod_sparse A;
		cholmod_common ch;
		cholmod_start(ch);
		A = cholmod_read_sparse(System.in, ch);
		if (A)
		{
			if (A.nrow != A.ncol || A.stype != 0
				||(!(A.xtype == CHOLMOD_REAL || A.xtype == CHOLMOD_COMPLEX)))
			{
				System.out.printf("invalid matrix\n");
			}
			else
			{
				klu_demo(A.nrow, A.p, A.i, A.x, A.xtype == CHOLMOD_REAL);
			}
			cholmod_free_sparse(A, ch);
		}
		cholmod_finish(ch);
		return(0);
	}

}
