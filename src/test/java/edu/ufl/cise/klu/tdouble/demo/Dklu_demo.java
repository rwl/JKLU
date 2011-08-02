package edu.ufl.cise.klu.tdouble.demo;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.mutable.MutableInt;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_version;
import edu.ufl.cise.klu.tdouble.Dklu_analyze;
import edu.ufl.cise.klu.tdouble.Dklu_defaults;
import edu.ufl.cise.klu.tdouble.Dklu_diagnostics;
import edu.ufl.cise.klu.tdouble.Dklu_factor;
import edu.ufl.cise.klu.tdouble.Dklu_free_numeric;
import edu.ufl.cise.klu.tdouble.Dklu_free_symbolic;
import edu.ufl.cise.klu.tdouble.Dklu_solve;
import edu.ufl.cise.klu.tdouble.io.MatrixInfo;
import edu.ufl.cise.klu.tdouble.io.MatrixSize;
import edu.ufl.cise.klu.tdouble.io.MatrixVectorReader;

/**
 * Read in a Matrix Market matrix(using CHOLMOD) and solve a linear system.
 */
public class Dklu_demo {

	private static void REAL (double[] X, int i, double v)
	{
		X[2*i] = v;
	}

	private static void IMAG (double[] X, int i, double v)
	{
		X[2*i + 1] = v;
	}

	private static double REAL (double[] X, int i)
	{
		return X[2*i] ;
	}

	private static double IMAG (double[] X, int i)
	{
		return X[2*i + 1] ;
	}

	private static double CABS (double[] X, int i)
	{
		return Math.sqrt(REAL(X, i) * REAL(X, i) + IMAG(X, i) * IMAG(X, i)) ;
	}

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
	 * @param lunz nnz(L+U+F)
	 * @param rnorm norm(b-A*x,1) / norm(A,1)
	 * @param Common default parameters and statistics
	 * @return 1 if successful, 0 otherwise
	 */
	public static int klu_backslash(int n, int[] Ap, int[] Ai, double[] Ax,
			int isreal, double[] B, double[] X, double[] R, MutableInt lunz,
			MutableDouble rnorm, KLU_common Common)
	{
		double anorm = 0, asum;
		KLU_symbolic Symbolic;
		KLU_numeric Numeric;
		int i, j, p;

		if (Ap == null || Ai == null || Ax == null || B == null || X == null ||
				B == null) return(0);

		/* ---------------------------------------------------------------------- */
		/* symbolic ordering and analysis */
		/* ---------------------------------------------------------------------- */

		Symbolic = Dklu_analyze.klu_analyze(n, Ap, Ai, Common);
		if (Symbolic == null) return(0);

		if (isreal != 0)
		{

			/* ------------------------------------------------------------------ */
			/* factorization */
			/* ------------------------------------------------------------------ */

			Numeric = Dklu_factor.klu_factor(Ap, Ai, Ax, Symbolic, Common);
			if (Numeric == null)
			{
				Dklu_free_symbolic.klu_free_symbolic(Symbolic, Common);
				return(0);
			}

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

			Dklu_diagnostics.klu_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common);
			Dklu_diagnostics.klu_condest(Ap, Ax, Symbolic, Numeric, Common);
			Dklu_diagnostics.klu_rcond(Symbolic, Numeric, Common);
			Dklu_diagnostics.klu_flops(Symbolic, Numeric, Common);
			lunz.setValue( Numeric.lnz + Numeric.unz - n +
				((Numeric.Offp != null) ? (Numeric.Offp [n]) : 0) );

			/* ------------------------------------------------------------------ */
			/* solve Ax=b */
			/* ------------------------------------------------------------------ */

			for(i = 0; i < n; i++)
			{
				X [i] = B [i];
			}
			Dklu_solve.klu_solve(Symbolic, Numeric, n, 1, X, Common);

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
			rnorm.setValue( 0 );
			for(i = 0; i < n; i++)
			{
				rnorm.setValue( MAX (rnorm.doubleValue(), Math.abs(R [i])) );
			}

			/* ------------------------------------------------------------------ */
			/* free numeric factorization */
			/* ------------------------------------------------------------------ */

			Dklu_free_numeric.klu_free_numeric(Numeric, Common);

		}
		else
		{
			throw new UnsupportedOperationException();

			/* ------------------------------------------------------------------ */
			/* statistics(not required to solve Ax=b) */
			/* ------------------------------------------------------------------ */

//			Numeric = DZklu_factor.klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
//			if (Numeric == null)
//			{
//				DZklu_free_symbolic.klu_free_symbolic(Symbolic, Common);
//				return(0);
//			}
//
//			/* ------------------------------------------------------------------ */
//			/* statistics */
//			/* ------------------------------------------------------------------ */
//
//			DZklu_diagnostics.klu_z_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common);
//			DZklu_diagnostics.klu_z_condest(Ap, Ax, Symbolic, Numeric, Common);
//			DZklu_diagnostics.klu_z_rcond(Symbolic, Numeric, Common);
//			DZklu_diagnostics.klu_z_flops(Symbolic, Numeric, Common);
//			lunz = Numeric.lnz + Numeric.unz - n +
//				((Numeric.Offp != null) ? (Numeric.Offp [n]) : 0);
//
//			/* ------------------------------------------------------------------ */
//			/* solve Ax=b */
//			/* ------------------------------------------------------------------ */
//
//			for(i = 0; i < 2*n; i++)
//			{
//				X [i] = B [i];
//			}
//			DZklu_solve.klu_z_solve(Symbolic, Numeric, n, 1, X, Common);
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
//			DZklu_free_numeric.klu_z_free_numeric(&Numeric, Common);
		}

		/* ---------------------------------------------------------------------- */
		/* free symbolic analysis, and residual */
		/* ---------------------------------------------------------------------- */

		Dklu_free_symbolic.klu_free_symbolic(Symbolic, Common);
		return(1);
	}

	/**
	 * Given a sparse matrix A, set up a right-hand-side and solve X = A\b.
	 */
	public static void klu_demo(int n, int[] Ap, int[] Ai, double[] Ax,
			int isreal)
	{
		MutableDouble rnorm = new MutableDouble();
		KLU_common Common = new KLU_common();
		double[] B, X, R;
		int i;
		MutableInt lunz = new MutableInt();

		System.out.printf("KLU: %s, version: %d.%d.%d\n",
				KLU_version.KLU_DATE, KLU_version.KLU_MAIN_VERSION,
				KLU_version.KLU_SUB_VERSION, KLU_version.KLU_SUBSUB_VERSION);

		/* ---------------------------------------------------------------------- */
		/* set defaults */
		/* ---------------------------------------------------------------------- */

		Dklu_defaults.klu_defaults(Common);

		/* ---------------------------------------------------------------------- */
		/* create a right-hand-side */
		/* ---------------------------------------------------------------------- */

		if (isreal != 0)
		{
			/* B = 1 +(1:n)/n */
//			B = klu_malloc(n, sizeof(Double), Common);
//			X = klu_malloc(n, sizeof(Double), Common);
//			R = klu_malloc(n, sizeof(Double), Common);
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
//			B = klu_malloc(n, 2 * sizeof(Double), Common);
//			X = klu_malloc(n, 2 * sizeof(Double), Common);
//			R = klu_malloc(n, 2 * sizeof(Double), Common);
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

		if (klu_backslash(n, Ap, Ai, Ax, isreal, B, X, R, lunz, rnorm, Common) == 0)
		{
			System.out.printf("KLU failed\n");
		}
		else
		{
			System.out.printf("n %d nnz(A) %d nnz(L+U+F) %d resid %g\n" +
				"recip growth %g condest %g rcond %g flops %g\n",
				n, Ap [n], lunz.intValue(), rnorm, Common.rgrowth, Common.condest,
				Common.rcond, Common.flops);
		}

		/* ---------------------------------------------------------------------- */
		/* free the problem */
		/* ---------------------------------------------------------------------- */

		if (isreal != 0)
		{
//			klu_free(B, n, sizeof(Double), Common);
//			klu_free(X, n, sizeof(Double), Common);
//			klu_free(R, n, sizeof(Double), Common);
			B = null;
			X = null;
			R = null;
		}
		else
		{
//			klu_free(B, 2*n, sizeof(Double), Common);
//			klu_free(X, 2*n, sizeof(Double), Common);
//			klu_free(R, 2*n, sizeof(Double), Common);
			B = null;
			X = null;
			R = null;
		}
		//System.out.printf("peak memory usage: %g bytes\n\n",(double)(Common.mempeak));
	}

	/**
	 * Read in a sparse matrix in Matrix Market format using CHOLMOD, and then
	 * solve Ax=b with KLU.  Note that CHOLMOD is only used to read the matrix.
	 */
	public static void main(String[] args) {

		FileReader fileReader;
		MatrixVectorReader reader;
		MatrixInfo info;
		MatrixSize size;

		try {
			fileReader = new FileReader(args[0]);
			reader = new MatrixVectorReader(fileReader);

			info = reader.readMatrixInfo();
			size = reader.readMatrixSize(info);

			if (size.numRows() != size.numColumns() || !info.isGeneral()
				||(!(info.isReal() || info.isComplex()))
				|| !info.isCoordinate())
			{
				System.out.printf("invalid matrix\n");
			}
			else
			{
				int nnz = size.numEntries();
				int[] i = new int[nnz];
				int[] p = new int[nnz];
				double[] x;

				if (info.isComplex()) {
					double[] xR = new double[nnz];
					double[] xI = new double[nnz];

					reader.readCoordinate(i, p, xR, xI);

					x = new double[2 * nnz];
					for (int j = 0; j < nnz; j++) {
						REAL(x, j, xR[j]);
						IMAG(x, j, xI[j]);
					}
				} else {
					x = new double[nnz];
					reader.readCoordinate(i, p, x);
				}


				klu_demo(size.numRows(), p, i, x, info.isReal() ? 1 : 0);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
