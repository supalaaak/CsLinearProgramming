using System;
using System.Data;
using System.Collections.Generic;
using System.IO;
using System.Linq;

public class CsLinearProgramming
{
	public static Double[][] Zdata;
	public static Double[] lowerBoundZero;
	public static string filePath;
	public static double kappa;
	public static double sigma;



	public class Bond
	{
		/*Global Variables*/
		public static double[] cashFlow;
		public static double[] cashDate;
		public static double static_spread;
		public static double[] callPrice;
		public static int numPeriod;

	}


	public static void setZeroMatrix(){

		StreamReader sr = new StreamReader(filePath);

		var lines = new List<double[]>();
		var lzero = new List<double> ();
		int Row = 0;
		Double num;
		while (!sr.EndOfStream)
		{

			String[] SLine = (sr.ReadLine().Split(','));
			if(double.TryParse(SLine[1], out num))
			{
				double[] Line =Array.ConvertAll<string,double>(SLine, double.Parse);
				lzero.Add (Line [0]);
				lines.Add(Line);
				Row++;
			}

		}
		Console.WriteLine (Row);
		Zdata = lines.ToArray();
		Console.WriteLine (Zdata.GetLength(0));

		lowerBoundZero = lzero.ToArray ();
	
	}
		
	public static double zero_map(double t){
		int idx = Array.IndexOf(lowerBoundZero,lowerBoundZero.Where(n => n <= t).Max());
		double[] izero = Zdata [idx];
		return izero [2] * t*t*t + izero [3] * t*t + izero [4] * t + izero [5];
	}

	public static double hull_white_r(double r0, double t1, double t2){

		int i1;
		double r1, r2,y,dy,d_r0,loga;
		r2 = zero_map (t2);
		double b = (1 - System.Math.Exp (-kappa * (t2 - t1))) / kappa;
		i1= Array.IndexOf(lowerBoundZero,lowerBoundZero.Where(n => n <= t1).Max());
		double[] izero = Zdata [i1];
		r1 = izero[2] * t1 * t1 * t1 +izero[3] * t1 * t1 + izero[4] * t1 + izero[5];
		y = 1 + r1 / 100;
		dy=(3 * izero[2] * t1 * t1 + 2 * izero[3] * t1 + izero[4]) / 100;
		int chk = Array.IndexOf(lowerBoundZero, t1);
		if (chk>= 0 && i1 + 1 < Zdata.GetLength (0)&&i1>0) {
			izero = Zdata [i1 + 1];
			dy=(dy+(3 * izero[2] * t1 * t1 + 2 * izero[3] * t1 + izero[4]) / 100)/2;
		}
		d_r0 = System.Math.Log (y) + t1 * dy / y;
		loga=t1*System.Math.Log(y)-t2*System.Math.Log(1+r2/100)+b*d_r0-sigma*sigma / 4 / kappa/kappa/kappa * System.Math.Pow(System.Math.Exp(-kappa * t2) - System.Math.Exp(-kappa * t1),2) * (System.Math.Exp(2 * kappa * t1) - 1);
		return System.Math.Exp(-b * r0 + loga);

	}

	public static double noAbritrage_bondPricing(int N,double[] cDate, double[] cFlow,double s_spread){

		double price = 0;
		for (int i = 0; i < N; i++) {
			double t = cDate [i];
			price = price + cFlow[i] / System.Math.Pow (1 + zero_map (t) / 100, t)*System.Math.Exp(-s_spread*t);
		}
		return price;
	}

	public static double HW_bondPricing(int N,double[] cDate, double[] cFlow,double s_spread,double t1,double r1){
	
		double price = 0;

		for (int i = 0; i < N; i++) {
			double t2 = cDate [i];
			if (t2 > t1) {
				price = price + cFlow [i]*hull_white_r(r1, t1, t2) * System.Math.Exp (-s_spread * (t2-t1));
			}
			}
		return price;
	}

	public static void  function1_fvec(double[] x, double[] fi, object obj)
	{
		//
		// this callback calculates
		// f0(x0,x1) = 100*(x0+3)^4,
		// f1(x0,x1) = (x1-3)^4
		//
		fi[0] = 10*System.Math.Pow(x[0]+3,2);
		fi[1] = System.Math.Pow(x[1]-3,2);
	}
	public static void function1_func(double[] x, ref double func, object obj)
	{
		// this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
		int fcall=Array.IndexOf(Bond.callPrice,Bond.callPrice.Where(n => n >0).First());
		func=System.Math.Pow(100- HW_bondPricing(Bond.numPeriod,Bond.cashDate, Bond.cashFlow,Bond.static_spread,Bond.cashDate[fcall],x[0]),2);
	}

	public static void function_cx_1_func(double[] c, double[] x, ref double func, object obj)
	{
		// this callback calculates f(c,x)=exp(-c0*sqr(x0))
		// where x is a position on X-axis and c is adjustable parameter
		int fcall=Array.IndexOf(Bond.callPrice,Bond.callPrice.Where(n => n >0).First());
		func=HW_bondPricing(Bond.numPeriod,Bond.cashDate, Bond.cashFlow,x[0],Bond.cashDate[fcall],c[0]);
	}
	public static void function2_cx_1_func(double[] c,double[] x, ref double func, object obj)
	{
		// this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
		func=100*System.Math.Pow(x[0]- CallaBond(c[0]),2);
	}

	public static double fitSpread_toCallPrice(double target_Price){

		double[,] x = new double[,]{{target_Price}};
		double[] y = new double[]{0};
		double[] c = new double[]{0.015};
		double[] bndl = new double[]{0.0};
		double[] bndu = new double[]{5.0};
		double epsf = 0;
		double epsx = 0.000001;
		int maxits = 0;
		int info;
		alglib.lsfitstate state;
		alglib.lsfitreport rep;
		double diffstep = 0.0001;

		alglib.lsfitcreatef(x, y, c, diffstep, out state);
		alglib.lsfitsetbc(state, bndl, bndu);
		alglib.lsfitsetcond(state, epsf, epsx, maxits);
		alglib.lsfitfit(state, function2_cx_1_func, null, null);
		alglib.lsfitresults(state, out info, out c, out rep);

		System.Console.WriteLine("Let's spread {0}", alglib.ap.format(c,5)); // EXPECTED: [1.0]

		return c [0];
	
	}

	public static double CallaBond(double s_spread){

		double straight_price = noAbritrage_bondPricing (Bond.numPeriod, Bond.cashDate, Bond.cashFlow, s_spread);
		int fcall=Array.IndexOf(Bond.callPrice,Bond.callPrice.Where(n => n >0).First());

		double[,] x = new double[,]{{s_spread}};
		double[] y = new double[]{Bond.callPrice[fcall]};
		double[] c = new double[]{0.015};
		double[] bndl = new double[]{0.0};
		double[] bndu = new double[]{5.0};
		double epsf = 0;
		double epsx = 0.000001;
		int maxits = 0;
		int info;
		alglib.lsfitstate state;
		alglib.lsfitreport rep;
		double diffstep = 0.0001;

		alglib.lsfitcreatef(x, y, c, diffstep, out state);
		alglib.lsfitsetbc(state, bndl, bndu);
		alglib.lsfitsetcond(state, epsf, epsx, maxits);
		alglib.lsfitfit(state, function_cx_1_func, null, null);
		alglib.lsfitresults(state, out info, out c, out rep);


		double t1 = Bond.cashDate [fcall];
		double sig_r0=sigma*System.Math.Sqrt((1-System.Math.Exp(-2*kappa*t1))/2/kappa);
		double z0 =System.Math.Pow(1+zero_map(t1)/100,-t1);
		double Option_price=0;
		double t2, K, z1, d1, d2,sb;
		double N = Bond.numPeriod;
		for (int i = fcall+1; i < N; i++) {
			t2 = Bond.cashDate [i];
			K = hull_white_r (c[0], t1, t2);
			sb=(1-System.Math.Exp(-kappa*(t2-t1)))/kappa*sig_r0;
			z1 = System.Math.Pow(1 + zero_map(t2)/100,-t2);
			d1 = (System.Math.Log (z1 / K / z0) + 0.5 * sb*sb) / sb ;
			d2 = d1 -  sb;
			Option_price = Option_price + Bond.cashFlow[i] * System.Math.Exp (-s_spread * t2) * (z1 * alglib.normaldistribution(d1) - K * z0 *alglib.normaldistribution (d2));

		}
		return  straight_price-Option_price;
	
	}
	public static int Main(string[] args)
	{
		//
		// This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
		//
		//     f0(x0,x1) = 10*(x0+3)^2
		//     f1(x0,x1) = (x1-3)^2
		//
		// using "V" mode of the Levenberg-Marquardt optimizer.
		//
		// Optimization algorithm uses:
		// * function vector f[] = {f1,f2}
		//
		// No other information (Jacobian, gradient, etc.) is needed.
		//
		kappa=1;
		sigma = 0.05;
		filePath="/Users/apple/projects/zero.csv";
		setZeroMatrix();
		Bond.cashFlow = new Double[]{ 5, 5, 5, 5, 5, 5, 5, 5, 5, 105 };
		Bond.cashDate = new Double[]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		Bond.callPrice = new Double[]{ 0, 0, 0, 0, 100, 100, 100,100, 100, 0 };
		Bond.numPeriod = 10;

		Bond.static_spread = 0.0271676179919139;
		double r0 = System.Math.Log (1 + zero_map (0) / 100);

		Console.WriteLine (noAbritrage_bondPricing(10,Bond.cashDate,Bond.cashFlow,Bond.static_spread));
		Console.WriteLine (HW_bondPricing(10,Bond.cashDate,Bond.cashFlow,Bond.static_spread,0,r0));

		double[] x = new double[]{0.015};
		double[] bndl = new double[]{0.0001};
		double[] bndu = new double[]{5};
		alglib.minbleicstate state;
		alglib.minbleicreport rep;

		//
		// These variables define stopping conditions for the optimizer.
		//
		// We use very simple condition - |g|<=epsg
		//
		double epsg = 0.000001;
		double epsf = 0;
		double epsx = 0;
		int maxits = 0;

		//
		// This variable contains differentiation step
		//
		double diffstep = 1.0e-6;

		//
		// Now we are ready to actually optimize something:
		// * first we create optimizer
		// * we add boundary constraints
		// * we tune stopping conditions
		// * and, finally, optimize and obtain results...
		//
		alglib.minbleiccreatef(x, diffstep, out state);
		alglib.minbleicsetbc(state, bndl, bndu);
		alglib.minbleicsetcond(state, epsg, epsf, epsx, maxits);
		alglib.minbleicoptimize(state, function1_func, null, null);
		alglib.minbleicresults(state, out x, out rep);


		System.Console.WriteLine("{0}", rep.terminationtype); // EXPECTED: 4
		System.Console.WriteLine("{0}", alglib.ap.format(x,5)); // EXPECTED: [-3,+3]
		Console.WriteLine (HW_bondPricing (10, Bond.cashDate, Bond.cashFlow, Bond.static_spread, 5, x[0]));
		int fcall=Array.IndexOf(Bond.callPrice,Bond.callPrice.Where(n => n >0).First());
		double t1 = Bond.cashDate [fcall];
		double sig_r0=sigma*System.Math.Sqrt((1-System.Math.Exp(-2*kappa*t1))/2/kappa);
		double z0 = hull_white_r (r0, 0, t1);//
		//z0=System.Math.Pow(1+zero_map(t1)/100,-t1);
		Double Option_price=0;
		double t2, K, z1, d1, d2,sb;
		for (int i = 5; i < 10; i++) {
			t2 = Bond.cashDate [i];
			K = hull_white_r (x[0], t1, t2);
			sb=(1-System.Math.Exp(-kappa*(t2-t1)))/kappa*sig_r0;
			z1 = System.Math.Pow(1 + zero_map(t2)/100,-t2);
			d1 = (System.Math.Log (z1 / K / z0) + 0.5 * sb*sb) / sb ;
			d2 = d1 - sb;
			Option_price = Option_price + Bond.cashFlow[i] * System.Math.Exp (-Bond.static_spread * t2) * (z1 * alglib.normaldistribution(d1) - K * z0 *alglib.normaldistribution (d2));
		}
		Console.WriteLine (Option_price);
		Console.WriteLine ("callable bond price {0}",CallaBond(0.0271676179919139));
		Console.WriteLine (CallaBond(fitSpread_toCallPrice(100)));

		System.Console.ReadLine();
		return 0;
	}
}