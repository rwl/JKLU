package edu.ufl.cise.klu.tdcomplex.util;

public class DoubleComplex {
	
	double[] component = new double[2];

	public DoubleComplex() {
	}

	public DoubleComplex(double[] component) {
		this.component = component;
	}
	
	public DoubleComplex(double re, double im) {
		this.component = new double[] {re, im};
	}
	
	double re() {
		return component[0];
	}
	
	double im() {
		return component[1];
	}

}
