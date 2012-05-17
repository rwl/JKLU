package edu.ufl.cise.klu.test;

import junit.framework.TestCase;
import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.tdcomplex.DZklu_common.DZklua;
import edu.ufl.cise.klu.tdcomplex.DZklu_internal;
import edu.ufl.cise.klu.tdcomplex.DZklu_version;

import static edu.ufl.cise.klu.tdcomplex.DZklu_defaults.klu_z_defaults;
import static edu.ufl.cise.klu.tdcomplex.DZklu_analyze.klu_z_analyze;
import static edu.ufl.cise.klu.tdcomplex.DZklu_factor.klu_z_factor;
import static edu.ufl.cise.klu.tdcomplex.DZklu_solve.klu_z_solve;

public class DZklu_simple extends TestCase {

	private static final double DELTA = 1e-09 ;

	private static int n = 5 ;
	private static int [ ] Ap = {0, 2, 5, 9, 10, 12} ;
	private static int [ ] Ai = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
	private static DZklua Ax = new DZklua(new double[] {2, 0, 3, 0, 3, 0, -1, 0, 4, 0, 4, 0, -3, 0, 1, 0, 2, 0, 2, 0, 6, 0, 1, 0}) ;
	private static DZklua b = new DZklua(new double[] {8, 0, 45, 0, -3, 0, 3, 0, 19, 0}) ;

	/**
	 * a simple KLU demo; solution is x = (1,2,3,4,5)
	 */
	public static void test_klu_simple() {
		int i;
		KLU_symbolic Symbolic;
		KLU_numeric Numeric;
		KLU_common Common = new KLU_common();

		DZklu_version.NPRINT = false ;
		DZklu_internal.NDEBUG = false ;

		klu_z_defaults (Common);
		Symbolic = klu_z_analyze (n, Ap, Ai, Common);
		Numeric = klu_z_factor (Ap, Ai, Ax, Symbolic, Common);
		klu_z_solve (Symbolic, Numeric, 5, 1, b, 0, Common);

		for (i = 0 ; i < n ; i++) {
			System.out.printf("x [%d] = %g\n", i, b.real(i)) ;
			assertEquals(i + 1.0, b.real(i), DELTA) ;
		}
	}

}
