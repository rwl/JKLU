package edu.ufl.cise.klu.tdcomplex;

public class DZklu_common {

	/**
	 *
	 * Complex array.
	 *
	 */
	public static class DZklua
	{

		/**
		 * numerical values
		 */
		public double[] x;

		public DZklua()
		{

		}

		public DZklua(double [] x)
		{
			this.x = x ;
		}

		/**
		 * Constructs an array of the given length.
		 */
		public DZklua(int len)
		{
			this.x = new double [2*len] ;
		}

		/**
		 *
		 * @param idx
		 * @return
		 */
		public double[] get(final int idx)
		{
			int offset = 2 * idx ;
			return new double [] {x [offset], x [offset + 1]} ;
		}

		public double real(final int idx) {
			return x [2 * idx] ;
		}

		public double imag(final int idx) {
			return x [(2 * idx) + 1] ;
		}

		public void real(final int idx, double re) {
			x [2 * idx] = re ;
		}

		public void imag(final int idx, double im) {
			x [(2 * idx) + 1] = im ;
		}

		/**
		 *
		 * @param idx
		 * @param val
		 */
		public void set(final int idx, final double [] val)
		{
			int offset = 2 * idx ;

			x [offset] = val [0] ;
			x [offset + 1] = val [1] ;
		}

		public void set(final int idx, final double re, final double im) {
			int offset = 2 * idx ;

			x [offset] = re ;
			x [offset + 1] = im ;
		}

		@Override
		public String toString() {
			String s = "DZcsa [" ;
			for (int i = 0; i < x.length; i+=2) {
				if (i != 0) s += ", " ;
				s += String.format("%g+j%g", x[i], x[i + 1]) ;
			}
			return s + "]" ;
		}
	}

}
