package edu.ufl.cise.klu.tdouble.demo;

import edu.ufl.cise.klu.tdouble.Dklu_solve;

public class Dklu_simple {

	/**
	 * a simple KLU demo; solution is x = (1,2,3,4,5)
	 */
	public static void main(String[] args) {
		int      n  = 5;
		int[]    Ap = {0, 2, 5, 9, 10, 12};
		int[]    Ai = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4};
		double[] Ax = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
		double[] b  = {8., 45., -3., 3., 19.};

	    Dklu_symbolic Symbolic;
	    Dklu_numeric Numeric;
	    Dklu_common Common;
	    int i;

	    Dklu_defaults.klu_defaults(Common);
	    Symbolic = Dklu_analyze.klu_analyze(n, Ap, Ai, Common);
	    Numeric = Dklu_factor.klu_factor(Ap, Ai, Ax, Symbolic, Common);
	    Dklu_solve.klu_solve(Symbolic, Numeric, 5, 1, b, Common);
	    Dklu_free_symbolic.klu_free_symbolic(Symbolic, Common);
	    Dklu_free_numeric.klu_free_numeric(Numeric, Common);
	    for (i = 0; i < n; i++)
	    	System.out.printf("x [%d] = %g\n", i, b [i]);
	}

}
