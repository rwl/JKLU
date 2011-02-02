package edu.ufl.cise.klu.tdcomplex.util;

public abstract class Util {
	
	public static final int INT_MAX = Integer.MAX_VALUE;
	
	public static final int KLU_OK = 0;
	public static final int KLU_SINGULAR = 1;            /* status > 0 is a warning, not an error */
	public static final int KLU_OUT_OF_MEMORY = -2;
	public static final int KLU_INVALID = -3;
	public static final int KLU_TOO_LARGE = -4;          /* integer overflow has occured */
	
	public static int MAX(int a, int b) {
		if (a > b) {
			return a;
		} else {
			return b;
		}
	}
	
	public static int MIN(int a, int b) {
		if (a < b) {
			return a;
		} else {
			return b;
		}
	}
	
	public static double MAX(double a, double b) {
		if (a > b) {
			return a;
		} else {
			return b;
		}
	}
	
	public static double MIN(double a, double b) {
		if (a < b) {
			return a;
		} else {
			return b;
		}
	}
	
	public static void PRINTF(String format,
			Object... args) {
		System.out.printf(format, args);
	}

	public static boolean SCALAR_IS_NAN(double x) {
		return Double.isNaN(x);
	}

	/* true if an integer (stored in double x) would overflow (or if x is NaN) */
	public static boolean INT_OVERFLOW(double x) {
		return (!((x) * (1.0+1e-8) <= (double) INT_MAX)) || SCALAR_IS_NAN (x);
	}
	
//	public static BYTES(type,n) (sizeof (type) * (n))
//	public static CEILING(b,u)  (((b)+(u)-1) / (u))
//	public static UNITS(type,n) (CEILING (BYTES (type,n), sizeof (Unit)))
//	public static DUNITS(type,n) (ceil (BYTES (type, (double) n) / sizeof (Unit)))
	
//	public static Object GET_POINTER(LU, Xip, Xlen, int Xi, Entry Xx, k, xlen) { 
//	    Unit xp = LU + Xip [k] ; 
//	    xlen = Xlen [k] ; 
//	    Xi = (int) xp ; 
//	    Xx = (Entry) (xp + UNITS (Int, xlen)) ; 
//	}

}
