// New parametrizations for sigma_diffractive
// by Appleby, Barlow, Molson, Serluca and Toader (ABMST).

//==========================================================================

// Necessary C++ libraries.
#include <cmath>
#include <iostream>
#include <iomanip>

// Simple shorthand.
using std::max;
using std::min;
using std::abs;
using std::cout;
using std::endl;
using std::fixed;
using std::scientific;
using std::setw;
using std::setprecision;

// Simple pow method.
inline double pow2(const double& x) {return x*x;}

//==========================================================================

// Single diffractive cross section according to
// Appleby, Barlow, Molson, Serluca and Toader (ABMST).

const double SPROTON = 0.8803544;
const double MPROTON = 0.9382720;
const double MPINEU  = 0.1349766;
const double SDAI[] = { 0.624529, 3.09088, 4.0, 177.217};
const double SDBI[] = { 2.5835, 4.51487, 3.03392, 5.86474};
const double SDCI[] = { 0.0, 0.186211, 10.0, 21.0029};
const double EPSI[] = { 0.08, -0.4525};
const double ALPP[] = { 0.25, 0.93};
const double POMS[] = { -0.25, -1.15, -0.05, 0.4, 0.5, 0.4597, 5.7575};
const double PIPS[] = { 13.63, 0.0808, 31.79, -0.4525, 14.4 };
const double MRES[] = { 1.44, 1.52, 1.68, 2.19};
const double WRES[] = { 0.325, 0.130, 0.140, 0.450};
const double CRES[] = { 3.07, 0.4149, 1.108, 0.9515};
const double FORM[] = { 2.79, 1., 0.71};

//--------------------------------------------------------------------------

// dsigma_sd/(dt dxi)(s, t, xi).
// iPart = -1: Merlin/Molson code as is.
//       = -2: Merlin/Molson, modified mCut according to article.
//       =  0: Torbjörn all.
//       =  1 - 7: Torbjörn separate components.

double dsigSDdtdxi( double s, double t, double xi, int iPart = 0) {

  // Code from Merlin for negative iPart.
  if (iPart < 0) {
        double x = xi;

	// Parameters for the resonance term in the background
	const double Mproton = 0.938272013;
	const double Mpion = 0.1349766;
	const double Mmin2 = pow(Mproton+Mpion,2);
	const double ml01 = 1.44;  //mass of the resonance P11 D13 G15 F17
	const double ml02 = 1.52;
	const double ml03 = 1.68;
	const double ml04 = 2.19;
	const double GammaL1 = 0.325; // width of the resonance
	const double GammaL2 = 0.13;
	const double GammaL3 = 0.14;
	const double GammaL4 = 0.45;
	const double cl01 = 3.07;   // coupling coefficient from the data fit
	const double cl02 = 0.4149;
	const double cl03 = 1.108  ;
	const double cl04 = 0.9515;
	//const double Mcut = 3; // this is chosen from the fit on the crosss section data
	//const double xi_c = pow(Mcut,2)/s;
	const double xi_th = Mmin2/s; // (M_p + M_pion)^2/s
	const double Mmin2bar = pow(Mproton-Mpion,2);
	if(x <= xi_th)
	{
		return 0;
	}

        // Modification of Merlin code: modified Mcut.
        double Mcut = 3.;
        if ( iPart == -2 && s > 4000.) Mcut = 3. + 0.6 * log( s / 4000.);
	double xi_c = pow(Mcut,2)/s;

	double cc[4][3];
	/*	cc[0][0] = 0.881148;	//A
		cc[0][1] = 3.94056;	//B
		cc[0][2] = 0.0220505;		//C

		//ppr
		cc[1][0] = 2.42997;
		cc[1][1] = 3.11514;
		cc[1][2] = 0.104746;

		//rrp
		cc[2][0] = 6.28648;
		cc[2][1] = 4.05376;
		cc[2][2] = 8.63609;

		//rrr
		cc[3][0] = 167.618;
		cc[3][1] = 11.5978;
		cc[3][2] = 54.256849;
	*/
	cc[0][0] = 0.624529;	//A
	cc[0][1] = 2.5835;	//B
	cc[0][2] = 0;		//C

	//ppr
	cc[1][0] = 3.09088;
	cc[1][1] = 4.51487;
	cc[1][2] = 0.186211;

	//rrp
	cc[2][0] = 4.0;
	cc[2][1] = 3.03392;
	cc[2][2] = 10.0;

	//rrr
	cc[3][0] = 177.217;
	cc[3][1] = 5.86474;
	cc[3][2] = 21.0029;



	const double q = sqrt((x*s-Mmin2)*(x*s-Mmin2bar)/(4*x*s));
	const double ql1 = sqrt( (pow(ml01,2) - Mmin2) * (pow(ml01,2) - Mmin2bar) /(4*pow(ml01,2)) );
	const double ql2 = sqrt( (pow(ml02,2) - Mmin2) * (pow(ml02,2) - Mmin2bar) /(4*pow(ml02,2)) );
	const double ql3 = sqrt( (pow(ml03,2) - Mmin2) * (pow(ml03,2) - Mmin2bar) /(4*pow(ml03,2)) );
	const double ql4 = sqrt( (pow(ml04,2) - Mmin2) * (pow(ml04,2) - Mmin2bar) /(4*pow(ml04,2)) );
	//std::cout << ql1 << '\t' << ql2 << '\t' << ql3 <<  '\t' << ql4 << std::endl;
	const double gammaL01 = GammaL1*pow(q/ql1,3)*((1 + 5*ql1)/(1 + 5*q));
	const double gammaL02 = GammaL2*pow(q/ql2,5)*pow(((1 + 5*ql2)/(1 + 5*q)),2);
	const double gammaL03 = GammaL3*pow(q/ql3,7)*pow(((1 + 5*ql3)/(1 + 5*q)),3);
	const double gammaL04 = GammaL4*pow(q/ql4,9)*pow(((1 + 5*ql4)/(1 + 5*q)),4);

	const double R = ( (cl01/x) * (ml01*gammaL01) / ( pow( (x*s - pow(ml01,2) ),2) + pow(ml01*gammaL01,2))
	                   +(cl02/x) * (ml02*gammaL02) / ( pow( (x*s - pow(ml02,2) ),2) + pow(ml02*gammaL02,2))
	                   +(cl03/x) * (ml03*gammaL03) / ( pow( (x*s - pow(ml03,2) ),2) + pow(ml03*gammaL03,2))
	                   +(cl04/x) * (ml04*gammaL04) / ( pow( (x*s - pow(ml04,2) ),2) + pow(ml04*gammaL04,2)) )
	                 *exp(13.5*(t + 0.05));// Normalization factors Sandy's note
	//* sqrt(565/s)*exp(13.5*(t + 0.05));// Normalization factors Sandy's note
//	double BRMatch = - 588.20982975 *exp(13.5*(t + 0.05))*(x - xi_th)/(xi_c - xi_th);
	const double BRMatch = -  ( (cl01/xi_c) * (ml01*gammaL01) / ( pow( (xi_c*s - pow(ml01,2) ),2) + pow(ml01*gammaL01,2))
	                            +(cl02/xi_c) * (ml02*gammaL02) / ( pow( (xi_c*s - pow(ml02,2) ),2) + pow(ml02*gammaL02,2))
	                            +(cl03/xi_c) * (ml03*gammaL03) / ( pow( (xi_c*s - pow(ml03,2) ),2) + pow(ml03*gammaL03,2))
	                            +(cl04/xi_c) * (ml04*gammaL04) / ( pow( (xi_c*s - pow(ml04,2) ),2) + pow(ml04*gammaL04,2)) )
	                       *exp(13.5*(t + 0.05)) *(x - xi_th)/(xi_c - xi_th);


//			   *sqrt(565/s)*exp(13.5*(t + 0.05))*(x - xi_th)/(xi_c - xi_th);


	//std::cout << "t = " << t << std::endl;
	if(t > -0.25)
	{
		//std::cout << "t less than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (0.4+0.5*t)*pow(s,0.08) * pow(xi_c,-1.08 -0.5*t)	//ppp
			                     +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			                     +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			                     +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			                     +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			                     * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor


			const double Aprimexi_c =(0.4+0.5*t)*pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t)	//ppp
			                         +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			                         +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			                         +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			                         +1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2)
			                         * ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))
			                            + (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor

			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);

			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;

		}
		else
		{
			return (0.4+0.5*t)*pow(s,0.08) * pow(x,-1.08 -0.5*t)	//ppp
			       +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			       +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			       +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			       +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			       * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}

	}
	else if(t > -1.15)
	{
		//std::cout << "t less than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(xi_c,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			                     +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			                     +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			                     +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			                     +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			                     * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor


			const double Aprimexi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			                          +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			                          +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			                          +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			                          +1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2)
			                          * ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))
			                             + (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor

			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);

			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;

		}
		else
		{
			return (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(x,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			       +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			       +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			       +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			       +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			       * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}

	}
	else
	{
		//std::cout << "t bigger or equal than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(xi_c,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			                     *(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			                     +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			                     +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			                     +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			                     +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			                     * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor


			const double Aprimexi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			                          *(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			                          +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			                          +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			                          +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			                          +1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2)
			                          * ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))
			                             + (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor

			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);

			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;

		}
		else
		{
			return (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(x,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			       *(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			       +(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			       +(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			       +(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			       +1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			       * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}

	}

  // End of code from Merlin. From now on Torbjörn's code.
  }

  // Vanishing cross section below the p + pi threshold.
  double spipSum = pow2( MPROTON + MPINEU);
  double spipDif = pow2( MPROTON - MPINEU);
  double m2X     = xi * s;
  if (m2X < spipSum) return 0.;

  // Calculate t range. Vanishing cross section if outside that range.
  double mu1     = SPROTON / s;
  double mu3     = xi;
  double rootv   = (1. - 4. * mu1) * (pow2(1. - mu1 - mu3) - 4. * mu1 * mu3);
  if (rootv <= 0.) return 0.;
  double tMin    = -0.5 * s * ( 1. - 3. * mu1 - mu3 + sqrt(rootv) );
  double tMax    = s * s * mu1 * pow2(mu3 - mu1) / tMin;
  if (t < tMin || t > tMax) return 0.;

  // Separation between low- and high-mass diffraction.
  // For low-mass evaluate P+R terms at the cut scale to allow matching.
  double m2Cut   = pow2( (s < 4000.) ? 3. : 3. + 0.6 * log( s / 4000.) );
  bool isHighM   = (m2X > m2Cut);
  double xiHigh  = (isHighM) ? xi : m2Cut / s;

  // Value of Pomeron/Reggeon trajectory.
  double alp0Pom = 1. + EPSI[0];
  double alp0Reg = 1. + EPSI[1];
  double alptPom = alp0Pom + ALPP[0] * t;
  double alptReg = alp0Reg + ALPP[1] * t;

  // PPP term, split by t range.
  /*
  double dsigPPP = t / (t + POMS[2])
    * pow( xiHigh, alp0Pom - 2. * alptPom) * pow( s, EPSI[0]);
  if (t > POMS[0]) dsigPPP *= POMS[3] + POMS[4] * t;
  else             dsigPPP *= (SDAI[0] * exp(SDBI[0] * t) + SDCI[0]);
  if (t < POMS[1]) dsigPPP *=  1. + POMS[5] * (POMS[1] - t)
    + POMS[6] * pow2(POMS[1] - t);
  */
  double dsigPPP = pow( xiHigh, alp0Pom - 2. * alptPom) * pow( s, EPSI[0]);
  if (t > POMS[0]) dsigPPP *= POMS[3] + POMS[4] * t;
  else dsigPPP *= ((SDAI[0] * exp(SDBI[0] * t) + SDCI[0])) * t / (t + POMS[2]);
  if (t < POMS[1]) dsigPPP *=  1. + POMS[5] * (POMS[1] - t)
    + POMS[6] * pow2(POMS[1] - t);

  // PPR, RRP and RRR terms.
  double dsigPPR = (SDAI[1] * exp(SDBI[1] * t) + SDCI[1])
    * pow( xiHigh, alp0Reg - 2. * alptPom) * pow( s, EPSI[1]);
  double dsigRRP = (SDAI[2] * exp(SDBI[2] * t) + SDCI[2])
    * pow( xiHigh, alp0Pom - 2. * alptReg) * pow( s, EPSI[0]);
  double dsigRRR = (SDAI[3] * exp(SDBI[3] * t) + SDCI[3])
    * pow( xiHigh, alp0Reg - 2. * alptReg) * pow( s, EPSI[1]);

  // Pion exchange term.
  double fForm   = (4. * SPROTON - FORM[0] * t)
                 / ( (4. * SPROTON - FORM[1] * t) * pow2(1. - t / FORM[2]) );
  double alptPi  = ALPP[1] * (t - pow2(MPINEU));
  double sigpip  = PIPS[0] * pow( xiHigh * s, PIPS[1])
                 + PIPS[2] * pow( xiHigh * s, PIPS[3]);
  double dsigPi  = PIPS[4] / (4. * M_PI) * (-t) / pow2(t - pow2(MPINEU))
                 * pow2(fForm) * pow( xiHigh, 1. - 2. * alptPi) * sigpip;

  // Done for high-mass diffraction.
  double dsigBkg = 0.;
  double dsigRes = 0.;
  double dsigTot = 0.;
  if (isHighM) dsigTot = dsigPPP + dsigPPR + dsigRRP + dsigRRR + dsigPi;

  // Low-mass diffraction: smoothly dampen high-mass contribution.
  else {
    double dsigCut = dsigPPP + dsigPPR + dsigRRP + dsigRRR + dsigPi;
    double derdCut = ( (alp0Pom - 2. * alptPom) * dsigPPP
      + (alp0Reg - 2. * alptPom) * dsigPPR
      + (alp0Pom - 2. * alptReg) * dsigRRP
      + (alp0Reg - 2. * alptReg) * dsigRRR
      + (1. - 2. * alptPi) * dsigPi ) / xiHigh;
    double xiThr   = spipSum / s;
    double derdDif = derdCut * (xiHigh - xiThr) - dsigCut;
    double coef2   = derdDif / pow2(xiHigh - xiThr);
    double coef1   = derdCut - 2. * derdDif / (xiHigh - xiThr);
    dsigBkg        = coef2 * pow2(xi - xiThr) + coef1 * (xi - xiThr);

    // Add low-mass resonances.
    double dsigSub = 0.;
    double qRef    = 0.5 * sqrt( (m2X - spipSum) * (m2X - spipDif) / m2X);
    for (int i = 0; i < 4; ++i) {
      double m2Now = pow2( MRES[i]);
      double qNow  = 0.5 * sqrt( (m2Now - spipSum) * (m2Now - spipDif) / m2Now);
      double mwNow = MRES[i] * WRES[i] * pow( qRef / qNow, 2 * i + 3)
                   * pow( (1. + 5. * qNow) / (1. + 5. * qRef), i + 1);
      dsigRes     += CRES[i] * mwNow / (pow2(m2X - m2Now) + pow2(mwNow));
      dsigSub     += CRES[i] * mwNow / (pow2(m2Cut - m2Now) + pow2(mwNow));
   }
    // Include Jacobian and t dependence, and subtract to vanish at cut.
    dsigSub       *= (xi - xiThr) / (xiHigh - xiThr);
    dsigRes        = (dsigRes / xi - dsigSub / xiHigh) * exp(13.5 * (t + 0.05));

    // Total low-mass contribution.
    dsigTot        = dsigBkg + dsigRes;
  }

  // Return sum, or for debug individual terms.
  if      (iPart == 0) return dsigTot;
  else if (iPart == 1 && isHighM) return dsigPPP;
  else if (iPart == 2 && isHighM) return dsigPPR;
  else if (iPart == 3 && isHighM) return dsigRRP;
  else if (iPart == 4 && isHighM) return dsigRRR;
  else if (iPart == 5 && isHighM) return dsigPi;
  else if (iPart == 6) return dsigBkg;
  else if (iPart == 7) return dsigRes;
  else                 return 0.;

}

//--------------------------------------------------------------------------

// Integrate over dt to obtain dsigma_sd/dxi(s, xi) in range tMin < t < tMax.

double dsigSDdxi( double s, double xi, double tMinIn = -1e10,
  double tMaxIn = 0., int iPart = 0) {

  // Number of integration points.
  int nPoints  = 200;

  // Calculate t range.
  double mu1   = SPROTON / s;
  double mu3   = xi;
  double rootv = (1. - 4. * mu1) * (pow2(1. - mu1 - mu3) - 4. * mu1 * mu3);
  if (rootv <= 0.) return 0.;
  double tMin  = -0.5 * s * ( 1. - 3. * mu1 - mu3 + sqrt(rootv) );
  double tMax  = s * s * mu1 * pow2(mu3 - mu1) / tMin;

  // Impose further t constraints from input. Check if range closed.
  tMin = max( tMin, tMinIn);
  tMax = min( tMax, tMaxIn);
  if (tMin > tMax) return 0.;

  // Prepare integration.
  double slope = -0.5 * log(xi);
  double etMin = exp(slope * tMin);
  double etMax = exp(slope * tMax);

  // Do integration by uniform steps in exp(slope * t).
  double dsigdxi= 0.;
  double y, yt, t;
  for (int i = 0; i < nPoints; ++i) {
    y          = (i + 0.5) / nPoints;
    yt         = etMin + y * (etMax - etMin);
    t          = log(yt) / slope;
    dsigdxi   += dsigSDdtdxi(s, t, xi, iPart) / yt;
  }

  // Normalize and done.
  dsigdxi *= (etMax - etMin) / (nPoints * slope);
  return dsigdxi;

}

//--------------------------------------------------------------------------

// Integrate over xi to obtain dsigma_sd/dt(s) in range xi_min < xi < xi_max.

double dsigSDdt( double s, double t, double xiMinIn = 0., double xiMaxIn = 1.,
  int iPart = 0) {

  // Density of integration points in ln(xi) for xi < xiDiv or else in xi.
  double xiDiv    = 0.1;
  double dlnxiRaw = 0.1;
  double dxiRaw   = 0.01;
  double mXmin    = MPROTON + MPINEU;

  // Approximate restrictions on xi range. Check it is not closed.
  double dsigdt   = 0.;
  double xiMin    = max( xiMinIn, mXmin * mXmin / s);
  double xiMax    = min( xiMaxIn, sqrt( -t / SPROTON) );
  if (xiMin >= xiMax) return 0.;

  // Integration in xi: check size of affected range and adjust dxi.
  if (xiMax > xiDiv) {
    double xiMinP = max( xiDiv, xiMin);
    int    nxi    = 2 + (xiMax - xiMinP) / dxiRaw;
    double dxi    = (xiMax - xiMinP) / nxi;
    for (int ixi = 0; ixi < nxi; ++ixi)
      dsigdt += dxi * dsigSDdtdxi( s, t, xiMinP + dxi * (ixi + 0.5), iPart);
  }

  // Integration in ln(xi): check size of affected range and adjust dlnxi.
  if (xiMin < xiDiv) {
    double xiMaxP = min( xiDiv, xiMax);
    int    nlnxi  = 2 + log( xiMaxP / xiMin) / dlnxiRaw;
    double dlnxi  = log( xiMaxP / xiMin) / nlnxi;
    for (int ilnxi = 0; ilnxi < nlnxi; ++ilnxi) {
      double xi   = xiMin * exp( dlnxi * (ilnxi + 0.5));
      dsigdt     += dlnxi * xi * dsigSDdtdxi( s, t, xi, iPart);
    }
  }

  // Done.
  return dsigdt;
}

//--------------------------------------------------------------------------

// Integrate over t and xi to obtain sigma_sd(s) in range  tMin < t < tMax
// and xi_min < xi < xi_max.

double sigSD( double s, double tMinIn = -1e10, double tMaxIn = 0.,
  double xiMinIn = 0., double xiMaxIn = 1., int iPart = 0) {

  // Density of integration points in ln(xi) for xi < xiDiv or else in xi.
  double xiDiv    = 0.1;
  double dlnxiRaw = 0.1;
  double dxiRaw   = 0.01;
  double mXmin    = 1.2;

  // Restrictions on xi range. Check it is not closed.
  double sig      = 0.;
  double xiMin    = max( xiMinIn, mXmin * mXmin / s);
  double xiMax    = min( xiMaxIn, sqrt( -tMinIn / SPROTON) );
  if (xiMin > xiMax) return 0.;

  // Integration in xi: check size of affected range and adjust dxi.
  if (xiMax > xiDiv) {
    double xiMinP = max( xiDiv, xiMin);
    int    nxi    = 2 + (xiMax - xiMinP) / dxiRaw;
    double dxi    = (xiMax - xiMinP) / nxi;
    for (int ixi = 0; ixi < nxi; ++ixi)
      sig        += dxi * dsigSDdxi( s, xiMinP + dxi * (ixi + 0.5),
                    tMinIn, tMaxIn, iPart);
  }

  // Integration in ln(xi): check size of affected range and adjust dlnxi.
  if (xiMin < xiDiv) {
    double xiMaxP = min( xiDiv, xiMax);
    int    nlnxi  = 2 + log( xiMaxP / xiMin) / dlnxiRaw;
    double dlnxi  = log( xiMaxP / xiMin) / nlnxi;
    for (int ilnxi = 0; ilnxi < nlnxi; ++ilnxi) {
      double xi   = xiMin * exp( dlnxi * (ilnxi + 0.5));
      sig        += dlnxi * xi * dsigSDdxi( s, xi, tMinIn, tMaxIn, iPart);
    }
  }

  // Done.
  return sig;
}

//--------------------------------------------------------------------------

// Alternative Integration over t and xi to obtain sigma_sd(s) in range
// tMin < t < tMax and xi_min < xi < xi_max.

double sigSDalt( double s, double tMin = -4., double tMax = 0.,
  double xiMin = 0., double xiMax = 1., int iPart = 0) {

  // Number of integration points in t and xi. Step size.
  int nt     = 500;
  int nxi    = 10000;
  double dt  = (tMax - tMin) / nt;
  double dxi = (xiMax - xiMin) / nxi;
  double sig = 0.;

  // Step through grid in t and xi.
  for (int it = 0; it < nt; ++it) {
    double t = tMin + it * dt;
    for (int ixi = 0; ixi < nxi; ++ixi) {
      double xi = xiMin + ixi * dxi;

      // Add up cross section.
      sig += dsigSDdtdxi( s, t, xi, iPart);
    }
  }

  // Normalize and done.
  sig *= dt * dxi;
  return sig;
}

//==========================================================================

int main() {

  // Table of integrated SD cross section.
  double ecmTab[14] = { 10., 14., 20., 30., 40., 60., 100., 200., 500.,
    1000., 2000., 5000., 10000., 20000. };
  double sigABMST[14] = { 0., 3.5, 4.5, 5.2, 5.6, 6.5, 6.8, 7.8, 9.3,
    0., 0., 0., 0., 0. };
  cout << "\n     eCM   Figure   Merlin   MerMod  Torbjorn  altInt " << endl;
  for (int iTab = 0; iTab < 14; ++iTab) {
    double sTab   = pow2(ecmTab[iTab]);
    double sigMer  = 2. * sigSD( sTab, -4.0, 0.0, 0.0,  0.05, -1);
    double sigMod  = 2. * sigSD( sTab, -4.0, 0.0, 0.0,  0.05, -2);
    double sigTor  = 2. * sigSD( sTab, -4.0, 0.0, 0.0,  0.05,  0);
    double sigTor2 = 2. * sigSDalt( sTab, -4.0, 0.0, 0.0,  0.05,  0);
    cout << fixed << setprecision(0) << setw(8) << ecmTab[iTab]
         << setprecision(3) << setw(9) << sigABMST[iTab] << setw(9)
         << sigMer << setw(9) << sigMod << setw(9) << sigTor
         << setw(9) << sigTor2 << endl;
  }

  return 0;
}
