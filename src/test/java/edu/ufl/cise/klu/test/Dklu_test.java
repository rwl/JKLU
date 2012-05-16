package edu.ufl.cise.klu.test;

import java.io.IOException;
import java.io.InputStream;

import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs;
import edu.emory.mathcs.csparsej.tdouble.test.Dcs_test;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_version;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_diagnostics.klu_rgrowth;
import static edu.ufl.cise.klu.tdouble.Dklu_diagnostics.klu_condest;
import static edu.ufl.cise.klu.tdouble.Dklu_diagnostics.klu_rcond;
import static edu.ufl.cise.klu.tdouble.Dklu_diagnostics.klu_flops;
import static edu.ufl.cise.klu.tdouble.Dklu_factor.klu_factor;
import static edu.ufl.cise.klu.tdouble.Dklu_solve.klu_solve;

/**
 * Read in a matrix and solve a linear system.
 */
public class Dklu_test extends Dcs_test {

	private static final String DIR = "matrix";

	private static final String ARROW = "arrow";

	private static final String IMPCOL_A = "impcol_a";

	private static final String WEST0156 = "west0156";

//	private static final String C1 = "1c";
//
//	private static final String ARROW_C = "arrowc";
//
//	private static final String C_TINA = "ctina";
//
//	private static final String GD99_CC = "GD99_cc";
//
//	private static final String ONE = "one";
//
//	private static final String ONE_C = "onec";
//
//	private static final String TWO = "two";
//
//	private static final String W156 = "w156";

	protected static InputStream get_stream(String name) {
		try
		{
			return Dklu_test.class.getResource(DIR + "/" + name).openStream() ;
		}
		catch (IOException e)
		{
			return (null) ;
		}
	}

	private static void REAL (double[] X, int i, double v)
	{
		X[2*i] = v;
	}

	private static void IMAG (double[] X, int i, double v)
	{
		X[2*i + 1] = v;
	}

//	private static double REAL (double[] X, int i)
//	{
//		return X[2*i] ;
//	}
//
//	private static double IMAG (double[] X, int i)
//	{
//		return X[2*i + 1] ;
//	}
//
//	private static double CABS (double[] X, int i)
//	{
//		return Math.sqrt(REAL(X, i) * REAL(X, i) + IMAG(X, i) * IMAG(X, i)) ;
//	}

	private static double MAX (double a, double b)
	{
		return (a) > (b) ? (a) : (b) ;
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
	 * @param lunz size 1, nnz(L+U+F)
	 * @param rnorm size 1, norm(b-A*x,1) / norm(A,1)
	 * @param Common default parameters and statistics
	 * @return 1 if successful, 0 otherwise
	 */
	public static int klu_backslash(int n, int[] Ap, int[] Ai, double[] Ax,
			boolean isreal, double[] B, double[] X, double[] R, int[] lunz,
			double[] rnorm, KLU_common Common)
	{
		double anorm = 0, asum;
		KLU_symbolic Symbolic;
		KLU_numeric Numeric;
		int i, j, p;

		if (Ap == null || Ai == null || Ax == null || B == null || X == null || B == null)
			return(0);

		/* ---------------------------------------------------------------------- */
		/* symbolic ordering and analysis */
		/* ---------------------------------------------------------------------- */

		Symbolic = klu_analyze (n, Ap, Ai, Common);
		if (Symbolic == null) return(0);

		if (isreal)
		{

			/* ------------------------------------------------------------------ */
			/* factorization */
			/* ------------------------------------------------------------------ */

			Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
			if (Numeric == null)
			{
				//klu_free_symbolic(Symbolic, Common);
				return(0);
			}

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

			klu_rgrowth (Ap, Ai, Ax, Symbolic, Numeric, Common);
			klu_condest (Ap, Ax, Symbolic, Numeric, Common);
			klu_rcond (Symbolic, Numeric, Common);
			klu_flops (Symbolic, Numeric, Common);
			lunz[0] = Numeric.lnz + Numeric.unz - n +
				(Numeric.Offp != null ? Numeric.Offp [n] : 0);

			/* ------------------------------------------------------------------ */
			/* solve Ax=b */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < n; i++)
			{
				X [i] = B [i];
			}
			klu_solve (Symbolic, Numeric, n, 1, X, 0, Common);

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
					asum += Math.abs(Ax [p]);
				}
				anorm = MAX (anorm, asum);
			}
			rnorm[0] = 0;
			for(i = 0; i < n; i++)
			{
				rnorm[0] = MAX (rnorm[0], Math.abs(R [i]));
			}

			/* ------------------------------------------------------------------ */
			/* free numeric factorization */
			/* ------------------------------------------------------------------ */

			//klu_free_numeric(Numeric, Common);
			Numeric = null;

		}
		else
		{
			throw new UnsupportedOperationException();

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

//			Numeric = klu_z_factor (Ap, Ai, Ax, Symbolic, Common);
//			if (Numeric == null)
//			{
//				klu_free_symbolic (Symbolic, Common);
//				return(0);
//			}
//
//			/* ------------------------------------------------------------------ */
//			/* statistics */
//			/* ------------------------------------------------------------------ */
//
//			klu_z_rgrowth (Ap, Ai, Ax, Symbolic, Numeric, Common);
//			klu_z_condest (Ap, Ax, Symbolic, Numeric, Common);
//			klu_z_rcond (Symbolic, Numeric, Common);
//			klu_z_flops (Symbolic, Numeric, Common);
//			lunz = Numeric.lnz + Numeric.unz - n +
//				(Numeric.Offp != null ? Numeric.Offp [n] : 0);
//
//			/* ------------------------------------------------------------------ */
//			/* solve Ax=b */
//			/* ------------------------------------------------------------------ */
//
//			for(i = 0; i < 2*n; i++)
//			{
//				X [i] = B [i];
//			}
//			klu_z_solve (Symbolic, Numeric, n, 1, X, Common);
//
//			/* ------------------------------------------------------------------ */
//			/* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
//			/* ------------------------------------------------------------------ */
//
//			for(i = 0; i < 2*n; i++)
//			{
//				R [i] = B [i];
//			}
//			for(j = 0; j < n; j++)
//			{
//				asum = 0;
//				for(p = Ap [j]; p < Ap [j+1]; p++)
//				{
//					/* R(i) -= A(i,j) * X(j) */
//					i = Ai [p];
//					REAL(R,i) -= REAL(Ax,p) * REAL(X,j) - IMAG(Ax,p) * IMAG(X,j);
//					IMAG(R,i) -= IMAG(Ax,p) * REAL(X,j) + REAL(Ax,p) * IMAG(X,j);
//					asum += CABS(Ax, p);
//				}
//				anorm = MAX(anorm, asum);
//			}
//			rnorm = 0;
//			for(i = 0; i < n; i++)
//			{
//				rnorm = MAX (rnorm, CABS(R, i));
//			}
//
//			/* ------------------------------------------------------------------ */
//			/* free numeric factorization */
//			/* ------------------------------------------------------------------ */
//
//			klu_z_free_numeric (&Numeric, Common);
		}

		/* ---------------------------------------------------------------------- */
		/* free symbolic analysis, and residual */
		/* ---------------------------------------------------------------------- */

		//klu_free_symbolic (Symbolic, Common);
		Symbolic = null;

		return (1);
	}

	/**
	 * Given a sparse matrix A, set up a right-hand-side and solve X = A\b.
	 */
	public static void klu_demo(int n, int[] Ap, int[] Ai, double[] Ax,
			boolean isreal, int[] lunz, double[] rnorm, KLU_common Common)
	{
		int i;
		double[] B, X, R;

		System.out.printf("KLU: %s, version: %d.%d.%d\n",
				KLU_version.KLU_DATE, KLU_version.KLU_MAIN_VERSION,
				KLU_version.KLU_SUB_VERSION, KLU_version.KLU_SUBSUB_VERSION);

		/* ---------------------------------------------------------------------- */
		/* set defaults */
		/* ---------------------------------------------------------------------- */

		klu_defaults (Common);

		/* ---------------------------------------------------------------------- */
		/* create a right-hand-side */
		/* ---------------------------------------------------------------------- */

		if (isreal)
		{
			/* B = 1 +(1:n)/n */
			B = new double[n];
			X = new double[n];
			R = new double[n];
			if (B != null)
			{
				for(i = 0; i < n; i++)
				{
					B [i] = 1 + ((double) i+1) / ((double) n);
				}
			}
		}
		else
		{
			/* real(B) = 1 +(1:n)/n, imag(B) = (n:-1:1)/n */
			B = new double[2 * n];
			X = new double[2 * n];
			R = new double[2 * n];
			if (B != null)
			{
				for(i = 0; i < n; i++)
				{
					REAL(B, i, 1 + ((double) i+1) / ((double) n));
					IMAG(B, i,     ((double) n-i) / ((double) n));
				}
			}
		}

		/* ---------------------------------------------------------------------- */
		/* X = A\b using KLU and print statistics */
		/* ---------------------------------------------------------------------- */

		if (klu_backslash (n, Ap, Ai, Ax, isreal, B, X, R, lunz, rnorm, Common) == 0)
		{
			System.out.printf("KLU failed\n");
		}
		else
		{
			System.out.printf("n %d nnz(A) %d nnz(L+U+F) %d resid %g\n" +
				"recip growth %g condest %g rcond %g flops %g\n",
				n, Ap [n], lunz[0], rnorm[0], Common.rgrowth, Common.condest,
				Common.rcond, Common.flops);
		}

		/* ---------------------------------------------------------------------- */
		/* free the problem */
		/* ---------------------------------------------------------------------- */

		B = null;
		X = null;
		R = null;

		System.out.printf("peak memory usage: %g bytes\n\n",(double)(Common.mempeak));
	}

	/**
	 * n 207 nnz(A) 572 nnz(L+U+F) 615 resid 6.98492e-10
	 * recip growth 0.00957447 condest 4.35093e+07 rcond 4.5277e-05 flops 259
	 */
	public void test_impcol_a() {
//		Dklu_version.NPRINT = false ;
//		Dklu_internal.NDEBUG = false ;

		KLU_common Common = new KLU_common();
		int[] lunz = new int[1];
		double[] rnorm = new double[1];

		InputStream in = get_stream (IMPCOL_A) ;
		Dproblem prob = get_problem (in, 0, 1) ;
		Dcs A = prob.A ;
		klu_demo (A.m, A.p, A.i, A.x, true, lunz, rnorm, Common) ;

		assertEquals(207, A.m) ;
		assertEquals(207, A.n) ;
		assertEquals(572, A.p[A.m]) ;
		assertEquals(615, lunz[0]) ;
		assertEquals(6.98492e-10, rnorm[0], 1e-14) ;
		assertEquals(0.00957447, Common.rgrowth, 1e-06) ;
		assertEquals(4.35093e+07, Common.condest, 1e+06) ;  // FIXME: improve assertion accuracy
		assertEquals(4.5277e-05, Common.rcond, 1e-08) ;
		assertEquals(259, Common.flops, 1e-03) ;
	}

	/**
	 * n 100 nnz(A) 298 nnz(L+U+F) 298 resid 1.77636e-15
	 * recip growth 0.0204082 condest 303 rcond 0.0204082 flops 297
	 */
	public void test_arrow() {
		KLU_common Common = new KLU_common() ;
		int[] lunz = new int[1] ;
		double[] rnorm = new double[1] ;

		InputStream in = get_stream (ARROW) ;
		Dproblem prob = get_problem (in, 0, 1) ;
		Dcs A = prob.A ;
		klu_demo (A.m, A.p, A.i, A.x, true, lunz, rnorm, Common) ;

		assertEquals(100, A.m) ;
		assertEquals(100, A.n) ;
		assertEquals(298, A.p[A.m]) ;
		assertEquals(298, lunz[0]) ;
		assertEquals(1.77636e-15, rnorm[0], 1e-18) ;
		assertEquals(0.0204082, Common.rgrowth, 1e-06) ;
		assertEquals(303, Common.condest, 1e-03) ;
		assertEquals(0.0204082, Common.rcond, 1e-06) ;
		assertEquals(297, Common.flops, 1e-03) ;
	}

	/**
	 * n 156 nnz(A) 371 nnz(L+U+F) 406 resid 1.04858e+06
	 * recip growth 0.0306751 condest 1.64225e+31 rcond 9.48528e-08 flops 188
	 */
	public void test_west0156() {
		KLU_common Common = new KLU_common();
		int[] lunz = new int[1];
		double[] rnorm = new double[1];

		InputStream in = get_stream (WEST0156) ;
		Dproblem prob = get_problem (in, 0, 1) ;
		Dcs A = prob.A ;
		klu_demo (A.m, A.p, A.i, A.x, true, lunz, rnorm, Common) ;

		assertEquals(156, A.m) ;
		assertEquals(156, A.n) ;
		assertEquals(371, A.p[A.m]) ;
		assertEquals(406, lunz[0]) ;
		assertEquals(1.04858e+06, rnorm[0], 1e+02) ;
		assertEquals(0.0306751, Common.rgrowth, 1e-06) ;
		assertEquals(1.64225e+31, Common.condest, 1e+26) ;
		assertEquals(9.48528e-08, Common.rcond, 1e-12) ;
		assertEquals(188, Common.flops, 1e-03) ;
	}

}
