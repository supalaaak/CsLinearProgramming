using System;
using System.Data;
using System.Collections.Generic;
using System.IO;
using System.Linq;

public class CsLinearProgramming
{
	/*Global HW Variables*/ 
	public static Double[][] Zdata;
	public static Double[] lowerBoundZero;
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


	public static void setZeroMatrix(String filePath){

		//
		// this callback calculates
		// set CB coefficient matrix and lowerBound
		// Zdata=[lowerbound upperbound coeffcient[1...4]
		//
		StreamReader sr = new StreamReader(filePath);

		var lines = new List<double[]>();
		var lzero = new List<double> ();
		Double num;
		while (!sr.EndOfStream)
		{

			String[] SLine = (sr.ReadLine().Split(','));
			if(double.TryParse(SLine[1], out num))
			{
				double[] Line =Array.ConvertAll<string,double>(SLine, double.Parse);
				lzero.Add (Line [0]);
				lines.Add(Line);
			}

		}
		Zdata = lines.ToArray();

		lowerBoundZero = lzero.ToArray ();
	
	}
		
	public static double zero_map(double t){
		//This callback calculates 
		//interest rates of the term t years based on coefficient matrix

		int idx = Array.IndexOf(lowerBoundZero,lowerBoundZero.Where(n => n <= t).Max());
		double[] izero = Zdata [idx];
		return izero [2] * t*t*t + izero [3] * t*t + izero [4] * t + izero [5];
	}

	public static double diff_zero_map(double t){
		
		//This callback calculates 
		//differentiate of the interest rates of the term t years based on coefficient matrix

		double dy=0;
		int idx = Array.IndexOf(lowerBoundZero,lowerBoundZero.Where(n => n <= t).Max());
		double[] izero = Zdata [idx];
		dy=(3*izero [2] * t*t + 2*izero [3] * t + izero [4] )/100;
		int chk = Array.IndexOf(lowerBoundZero, t);
		if (chk> -1 && idx + 1 < Zdata.GetLength (0) && idx>0) {
			izero = Zdata [idx + 1];
			dy=(dy+(3 * izero[2] * t * t + 2 * izero[3] * t + izero[4]) / 100)/2;
		}
		return dy;
	}



	public static double hull_white_r(double r0, double t1, double t2){

		//This callback calculates the present value of a dollar to be paid in t2 with initial time t1
		//where r0 is initial interest rate with respect to t1
		//function hull_white_r=exp(-b*r0+loga)
		//  b=(1-exp(-kappa*(t2-t1))/kappa;
		//  loga=log(Z(0,t2)/Z(0,t1))-b*d/dt(log(Z(0,t1))-sig^2/4/kappa^3*(exp(-kappa*t2)
		//       -exp(-kappa*t1))^2*(exp(2*kappa*t1)-1);

		double r1, r2,y,dy,d_r0,loga;
		r1 = zero_map (t1);
		r2 = zero_map (t2);
		double b = (1 - System.Math.Exp (-kappa * (t2 - t1))) / kappa;
		y = 1 + r1 / 100;
		dy = diff_zero_map (t1);
		d_r0 = System.Math.Log (y) + t1 * dy / y;
		loga=t1*System.Math.Log(y)-t2*System.Math.Log(1+r2/100)+b*d_r0-sigma*sigma / 4 / kappa/kappa/kappa * System.Math.Pow(System.Math.Exp(-kappa * t2) - System.Math.Exp(-kappa * t1),2) * (System.Math.Exp(2 * kappa * t1) - 1);
		return System.Math.Exp(-b * r0 + loga);

	}

	public static double noAbitrage_bondPricing(int N,double[] cDate, double[] cFlow,double s_spread){

		/*
		This callback calculates price of bond using no Arbitrage model
		function noAbritprice=sum c(i)/(1+r(ti)^ti*exp(-s_spread*ti)
		where N is the number of cashflow, 
		cDate is effective periodic time,
		cFlow is effective cash flow, 
		s_spread is static spread
        */

		double price = 0;
		for (int i = 0; i < N; i++) {
			double t = cDate [i];
			price = price + cFlow[i] / System.Math.Pow (1 + zero_map (t) / 100, t)*System.Math.Exp(-s_spread*t);
		}
		return price;
	}

	public static double HW_bondPricing(int N,double[] cDate, double[] cFlow,double s_spread,double t1,double r1){
	
		//This callback calculates price of bond using HW model from present time t1
		//function noAbritprice=sum c(i)*hull_white_r(r1,t1,ti)*exp(-spread*ti)
		// where N is the number of cashflow, cDate is effective periodic time,
		// cFlow is effective cash flow, s_spread is static spread
		//  and r1 is perspective initial interest rate

		double price = 0;

		for (int i = 0; i < N; i++) {
			double t2 = cDate [i];
			if (t2 > t1) {
				price = price + cFlow [i]*hull_white_r(r1, t1, t2) * System.Math.Exp (-s_spread * (t2-t1));
			}
			}
		return price;
	}



	public static void function_cx_1_func(double[] c, double[] x, ref double func, object obj)
	{
		// this function is used for optimization to solve for interest rate c[0] at first call date
		// ....that makes sense for the bond to be early redeemed.
		// this callback calculates bond price by using HW model 
		// where x[0] is static spread and c[0] is  interest rate with respect first time call

		int fcall=Array.IndexOf(Bond.callPrice,Bond.callPrice.Where(n => n >0).First());
		func=HW_bondPricing(Bond.numPeriod,Bond.cashDate, Bond.cashFlow,x[0],Bond.cashDate[fcall],c[0]);
	}
	public static void function2_cx_1_func(double[] c,double[] x, ref double func, object obj)
	{
		// this function is used for optimization to solve for static spread c[0] that makes the callable bond price
		// equal to the target price x[0]
		//  obj  min target price -callable bond price
		//   solve for c[0];

 		func=100*System.Math.Pow(x[0]- CallaBond(c[0]),2);
	}

	public static double fitSpread_toCallPrice(double target_Price){

		// this callback returns static_spread that makes callable price equal to target_price
		// by using Nonlinear least squares fitting of alglib package
		/*
		  Int this function we demonstrate fitting by
		    target price x[0]= Callabond(c[0])
		  subject to bound constraints
		     0<c[0]<=5
		  Our variables are:
           c[0]=static spread
		   x[0]=target price;

        */
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


		return c [0];
	
	}
	public static double CallaBond(){

		// this callback returns callable bond price
		// by using Nonlinear least squares fitting of alglib package and
		// alglib.normaldistribution for cumulative distribution function (cdf) of normal
		// where s_spread is static spread
		double s_spread=Bond.static_spread;
		double straight_price = noAbitrage_bondPricing (Bond.numPeriod, Bond.cashDate, Bond.cashFlow, s_spread);
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

	public static double CallaBond(double s_spread){
		
		// this callback returns callable bond price
		// by using Nonlinear least squares fitting of alglib package and
		// alglib.normaldistribution for cumulative distribution function (cdf) of normal
	    // where s_spread is static spread
      
		double straight_price = noAbitrage_bondPricing (Bond.numPeriod, Bond.cashDate, Bond.cashFlow, s_spread);
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
		//initialize HW parameters kappa,sigma, CB coefficients
		kappa=1;
		sigma = 0.05;
		String filePath=@"/Users/apple/projects/zero.csv";
		setZeroMatrix(filePath);
		//Intialized r0 (initial interest rate at time 0
		double r0 = System.Math.Log (1 + zero_map (0) / 100);


		//Initialize Bond cashFlow, cashDate,call_Price,num_periods,static spread
		Bond.cashFlow = new Double[]{ 5, 5, 5, 5, 5, 5, 5, 5, 5, 105 };
		Bond.cashDate = new Double[]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		Bond.callPrice = new Double[]{ 0, 0, 0, 0, 100, 100, 100,100, 100, 100 };
		Bond.numPeriod = 10;
		Bond.static_spread = 0.02;

		Console.WriteLine ("With statc spread 0.027 the straight bond price is {0}",noAbitrage_bondPricing(10,Bond.cashDate,Bond.cashFlow,Bond.static_spread));
		Console.WriteLine ("With statc spread 0.027 the straight bond price using HW Model is {0}",HW_bondPricing(10,Bond.cashDate,Bond.cashFlow,Bond.static_spread,0,r0));
		Console.WriteLine ("and Callable bond price is {0}",CallaBond());
		Console.WriteLine ("If the market price of a callable bond is 100 the static spread must be {0}",fitSpread_toCallPrice(100));

		return 0;
	}
}
