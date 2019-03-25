// New parametrizations for sigma_tot and dsigma_el/dt
// by Appleby, Barlow, Molson, Serluca and Toader (ABMST),
// and by the Review of Particle Physics (RPP).

//==========================================================================

// Code normally found in PythiaStdlib.h and PythiaComplex.h.

/*
#include <cmath>
#include <complex>
typedef std::complex<double> complex;
#include <iostream>
#include <iomanip>
using std::abs;
using std::cout;
using std::endl;
using std::fixed;
using std::scientific;
using std::setw;
using std::setprecision;
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
*/

#include "Pythia8/Basics.h"
#include "Pythia8/PythiaComplex.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Complex Bessel functions J0 and J1.

complex besJ0( complex x) {
  int mMax    = 5. + 5. * abs(x);
  complex z   = 0.25 * x * x;
  complex term = 1.;
  complex sum  = term;
  for (int m = 1; m < mMax; ++m) {
    term *= - z / double(m * m);
    sum  += term;
  }
  return sum;
}
complex besJ1( complex x) {
  int mMax    = 5. + 5. * abs(x);
  complex z   = 0.25 * x * x;
  complex term = 0.5 * x;
  complex sum  = term;
  for (int m = 1; m < mMax; ++m) {
    term *= - z / double(m * (m+1));
    sum  += term;
  }
  return sum;
}

//--------------------------------------------------------------------------

// Common values.
const double SPROTON = 0.8803544;
const double HBARC2  = 0.38938;

// Common method.
complex sModAlp( double sMod, double alpha) {
  return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha);
}

//==========================================================================

// Total and elastic cross section according to
// Appleby, Barlow, Molson, Serluca and Toader (ABMST).

const double EPSI[]  = { 0.106231, 0.0972043, -0.510662, -0.302082};
const double ALPP[]  = { 0.0449029, 0.278037, 0.821595, 0.904556};
const double NORM[]  = { 228.359, 193.811, 518.686, 10.7843};
const double SLOPE[] = { 8.38, 3.78, 1.36};
const double FRACS[] = { 0.26, 0.56, 0.18};
const double TRIG[]  = { 0.3, 5.03};
const double LAMBDA  = 0.521223;

//--------------------------------------------------------------------------

// Amplitude.

complex amplitudeABMST( double s, double t, double ispp = true) {

  // Common values.
  double snu  = s - 2. * SPROTON + 0.5 * t;
  double ampt = FRACS[0] * exp(SLOPE[0] * t) + FRACS[1] * exp(SLOPE[1] * t)
              + FRACS[2] * exp(SLOPE[2] * t);
  complex amp[6], l2p[4], ll2p[4], d2p[4][3];

  // Two Pomeron and even and odd Reggeon exchange.
  for (int i = 0; i < 4; ++i)
    amp[i] = ((i < 3) ? complex(-NORM[i], 0.) : complex( 0., NORM[i]))
           * ampt * sModAlp( ALPP[i] * snu, 1. + EPSI[i] + ALPP[i] * t);

  // Two-pomeron exchange.
  amp[4] = complex(0., 0.);
  for (int i = 0; i < 4; ++i) {
    l2p[i]  = ALPP[i] * complex( log(ALPP[i] * snu), -0.5 * M_PI);
    ll2p[i] = (1. + EPSI[i]) * l2p[i] / ALPP[i];
    for (int k = 0; k < 3; ++k) d2p[i][k] = SLOPE[k] + l2p[i];
  }
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
  for (int k = 0; k < 3; ++k)
  for (int l = 0; l < 3; ++l) {
    complex part = NORM[i] * NORM[j] * exp( ll2p[i] + ll2p[j] )
                 * exp( t * d2p[i][k] * d2p[j][l] / (d2p[i][k] + d2p[j][l]) )
                 * FRACS[k] * FRACS[l] / (d2p[i][k] + d2p[j][l]);
    if (i == 3) part *= complex( 0., 1.);
    if (j == 3) part *= complex( 0., 1.);
    amp[4]      += part;
  }
  amp[4]        *= LAMBDA * complex( 0., 1.) / (16. * M_PI * snu);

  // Triple-gluon exchange.
  amp[5] = sqrt(16. * M_PI / HBARC2) * TRIG[0] * ((t < -TRIG[1])
         ? 1. / pow4(t) :  exp(4. + 4. * t / TRIG[1]) / pow4(TRIG[1]));

  // Add up contributions.
  complex ampSum = (amp[0] + amp[1] + amp[2] + ((ispp) ? -amp[3] : amp[3])
    + amp[4]) / snu + ((ispp) ? amp[5] : -amp[5]);
  return ampSum;

}

//--------------------------------------------------------------------------

// Total cross section.

double sigmaABMST(double s, double ispp = true) {
  return HBARC2 * imag(amplitudeABMST( s, 0., ispp));
}

//--------------------------------------------------------------------------

// The rho parameter.

double rhoABMST(double s, double ispp = true) {
  complex amp = amplitudeABMST( s, 0., ispp);
  return real(amp) / imag(amp);
}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double dsigmadtABMST(double s, double t, double ispp = true) {
  return HBARC2 * pow2(abs(amplitudeABMST( s, t, ispp))) / (16. * M_PI);
}

//--------------------------------------------------------------------------

// Slope parameter , defined between t = 0 and t = -0.01 by default.

double bSlopeABMST(double s, double ispp = true, double t2 = -0.01) {
  return log( dsigmadtABMST(s, t2, ispp) / dsigmadtABMST(s, 0., ispp) ) / t2;
}

//--------------------------------------------------------------------------

// Total elastic cross section, by integration or optical theorem.

double sigmaelABMST(double s, double ispp = true, bool integrate = true) {
  double sigel = 0.;
  if (integrate) {
    int nPoints = 1000;
    double samplingSlope = 10.;
    for (int i = 0; i < nPoints; ++i) {
      double y = (i + 0.5) / nPoints;
      double t = log(y) / samplingSlope;
      sigel += dsigmadtABMST(s, t, ispp) / y;
    }
    sigel /= (nPoints * samplingSlope);
  } else {
    sigel = pow2(sigmaABMST(s, ispp)) * (1. + pow2(rhoABMST(s, ispp)))
      / (HBARC2 * 16. * M_PI * bSlopeABMST(s, ispp));
  }
  return sigel;
}

//====================================================================

// Total and elastic cross section according to 2014 RPP.

const double EPS1RPP[]  = { 1., 0.614, 0.444, 1., 1., 1.};
const double ALPPRPP[]  = { 0.151, 0.8, 0.8, 0.947};
const double NORMRPP[] = { 0.2478, 0.0078, 11.22, -0.150, 148.4, -26.6,
  -1.5, -0.0441, 0., 0.686, -3.82, -8.60, 64.1, 99.1, -58.0, 9.5};
const double BRPP[] = { 3.592, 0.622, 5.44, 0.205, 5.643, 1.92, 0.41, 0.,
  0., 3.013, 2.572, 12.25, 2.611, 11.28, 1.27};
const double KRPP[] = { 0.3076, 0.0998, 1.678, 0.190, -26.1};


//--------------------------------------------------------------------------

// Amplitude.

complex amplitudeRPP(double s, double t, double ispp = true) {

  // Modified s-related values.
  double  shat   = s - 2. * SPROTON + 0.5 * t;
  complex stlog  = complex( log(shat), -0.5 * M_PI);
  complex taut   = sqrt(abs(t)) * stlog;

  // Trajectories.
  double aP      = EPS1RPP[0] + ALPPRPP[0] * t;
  double aRpos   = EPS1RPP[1] + ALPPRPP[1] * t;
  double aRneg   = EPS1RPP[2] + ALPPRPP[2] * t;
  double aO      = EPS1RPP[3] + ALPPRPP[3] * t;
  double aOP     = EPS1RPP[4] + ALPPRPP[0] * ALPPRPP[3] * t
                 / (ALPPRPP[0] + ALPPRPP[3]);
  double aPP     = EPS1RPP[5] + 0.5 * ALPPRPP[0] * t;
  double aRPpos  = EPS1RPP[1] + ALPPRPP[0] * ALPPRPP[1] * t
                 / (ALPPRPP[0] + ALPPRPP[1]);
  double aRPneg  = EPS1RPP[2] + ALPPRPP[0] * ALPPRPP[2] * t
                 / (ALPPRPP[0] + ALPPRPP[2]);

  // Even terms.
  complex besArg = KRPP[0] * taut;
  complex besJ0n = besJ0(besArg);
  complex besJ1n = besJ1(besArg);
  complex besRat = (abs(besArg) < 0.01) ? 1. : 2. * besJ1n / besArg;
  complex fPosH  = complex( 0., shat) * (NORMRPP[0] * besRat
                 * exp(BRPP[0] * t) * stlog * stlog
                 + NORMRPP[1] * besJ0n * exp(BRPP[1] * t) * stlog
                 + NORMRPP[2] * (besJ0n - besArg * besJ1n) * exp(BRPP[2] * t));
  complex fPosP  = -NORMRPP[3] * exp(BRPP[3] * t) * sModAlp( shat, aP);
  complex fPosPP = (-NORMRPP[4] / stlog) * exp(BRPP[4] * t)
                 * sModAlp( shat, aPP);
  complex fPosR  = -NORMRPP[5] * exp(BRPP[5] * t) * sModAlp( shat, aRpos);
  complex fPosRP = t * (NORMRPP[6] / stlog) * exp(BRPP[6] * t)
                 * sModAlp( shat, aRPpos);
  complex nPos   = complex( 0., -shat) * NORMRPP[7] * stlog * t
                 * pow( 1. - t / KRPP[2], -5.);
  complex fPos   = fPosH + fPosP + fPosPP + fPosR + fPosRP + nPos;

  // Odd terms.
  complex fNegMO = shat * (NORMRPP[9] * cos(KRPP[1] * taut) * exp(BRPP[9] * t)
                 * stlog + NORMRPP[10] * exp(BRPP[10] * t));
  complex fNegO  = complex( 0., NORMRPP[11]) * exp(BRPP[11] * t)
                 * sModAlp( shat, aO) * (1. + KRPP[4] * t);
  complex fNegOP = (complex( 0., NORMRPP[12]) / stlog) * exp(BRPP[12] * t)
                 * sModAlp( shat, aOP);
  complex fNegR  = complex( 0., -NORMRPP[13]) * exp(BRPP[13] * t)
                 * sModAlp( shat, aRneg);
  complex fNegRP = t * (complex( 0., -NORMRPP[14]) / stlog) * exp(BRPP[14] * t)
                 * sModAlp( shat, aRPneg);
  complex nNeg   = -shat * NORMRPP[15] * stlog * t
                 * pow( 1. - t / KRPP[3], -5.);
  complex fNeg   = fNegMO + fNegO + fNegOP + fNegR + fNegRP + nNeg;

  // Combine and return.
  complex ampSum = ispp ? fPos + fNeg : fPos - fNeg;
  return  ampSum;

}

//--------------------------------------------------------------------------

// Total cross section.

double sigmaRPP(double s, double ispp = true) {
  return imag(amplitudeRPP( s, 0., ispp)) / sqrt(s * ( s - 4. * SPROTON));
}

//--------------------------------------------------------------------------

// The rho parameter.

double rhoRPP(double s, double ispp = true) {
  complex amp = amplitudeRPP( s, 0., ispp);
  return real(amp) / imag(amp);
}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double dsigmadtRPP(double s, double t, double ispp = true) {
  return pow2(abs(amplitudeRPP( s, t, ispp)))
         / (16. * M_PI * HBARC2 * s * ( s - 4. * SPROTON));
}

//--------------------------------------------------------------------------

// Slope parameter , defined between t = 0 and t = 0.01 by default.

double bSlopeRPP(double s, double ispp = true, double t2 = -0.01) {
  return log( dsigmadtRPP(s, t2, ispp) / dsigmadtRPP(s, 0., ispp) ) / t2;
}

//--------------------------------------------------------------------------

// Total elastic cross section, by integration or optical theorem.

double sigmaelRPP(double s, double ispp = true, bool integrate = true) {
  double sigel = 0.;
  if (integrate) {
    int nPoints = 1000;
    double samplingSlope = 10.;
    for (int i = 0; i < nPoints; ++i) {
      double y = (i + 0.5) / nPoints;
      double t = log(y) / samplingSlope;
      sigel += dsigmadtRPP(s, t, ispp) / y;
    }
    sigel /= (nPoints * samplingSlope);
  } else {
    sigel = pow2(sigmaRPP(s, ispp)) * (1. + pow2(rhoRPP(s, ispp)))
      / (HBARC2 * 16. * M_PI * bSlopeRPP(s, ispp));
  }
  return sigel;
}

//==========================================================================

int main() {

  // Points in energy for sigma_tot and in t for dsigma_el/dt.
  //double eCMval[14] = { 5., 10., 20., 50., 100., 200., 500., 1000.,
  //  2000., 5000., 7000., 8000., 10000., 13000.};
  double eCMval[4] = { 10., 100., 1000., 10000.};
  int nE = 4;
  //double tval[18] = { 0., -0.0005, -0.001, -0.002, -0.005, -0.01, -0.02,
  //  -0.05, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -1., -2., -5.};
  //int nT = 18;

  // Energy dependence of sigma_tot.
  cout << "\n   eCM        sigma_pp           sigma_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM        = eCMval[iE];
    double s          = eCM * eCM;
    double sigPPABMST = sigmaABMST(s, true);
    double sigPbABMST = sigmaABMST(s, false);
    double sigPPRPP   = sigmaRPP(s, true);
    double sigPbRPP   = sigmaRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << sigPPABMST << setw(10) << sigPPRPP
         << setw(10) << sigPbABMST << setw(10) << sigPbRPP << endl;
  }

  // Energy dependence of the rho parameter.
  cout << "\n   eCM         rho_pp             rho_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM        = eCMval[iE];
    double s          = eCM * eCM;
    double rhoPPABMST = rhoABMST(s, true);
    double rhoPbABMST = rhoABMST(s, false);
    double rhoPPRPP   = rhoRPP(s, true);
    double rhoPbRPP   = rhoRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << rhoPPABMST << setw(10) << rhoPPRPP
         << setw(10) << rhoPbABMST << setw(10) << rhoPbRPP << endl;
  }

  // Energy dependence of sigma_elastic; nmerical integration.
  cout << "\n   sigma_elastic by numerical integration"
       << "\n   eCM       sigmael_pp         sigmael_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM        = eCMval[iE];
    double s          = eCM * eCM;
    double sigPPABMST = sigmaelABMST(s, true, true);
    double sigPbABMST = sigmaelABMST(s, false, true);
    double sigPPRPP   = sigmaelRPP(s, true, true);
    double sigPbRPP   = sigmaelRPP(s, false, true);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << sigPPABMST << setw(10) << sigPPRPP
         << setw(10) << sigPbABMST << setw(10) << sigPbRPP << endl;
  }

  // Energy dependence of sigma_elastic; optical theorem + fix slope.
  cout << "\n   sigma_elastic via optical theorem + fix slope"
       << "\n   eCM       sigmael_pp         sigmael_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM        = eCMval[iE];
    double s          = eCM * eCM;
    double sigPPABMST = sigmaelABMST(s, true, false);
    double sigPbABMST = sigmaelABMST(s, false, false);
    double sigPPRPP   = sigmaelRPP(s, true, false);
    double sigPbRPP   = sigmaelRPP(s, false, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << sigPPABMST << setw(10) << sigPPRPP
         << setw(10) << sigPbABMST << setw(10) << sigPbRPP << endl;
  }

  // Approximate t slopes.
  cout << "\n   eCM          B_pp             B_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double bPPABMST = bSlopeABMST(s, true);
    double bPbABMST = bSlopeABMST(s, false);
    double bPPRPP   = bSlopeRPP(s, true);
    double bPbRPP   = bSlopeRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << bPPABMST << setw(10) << bPPRPP
         << setw(10) << bPbABMST << setw(10) << bPbRPP << endl;
  }

  /*
  // Energy dependence of dsigma/dt at some fixed t; here -0.1.
  double tNow = -0.1;
  cout << "\n   dsigma/dt at t = " << fixed << setprecision(4) << tNow
       << "\n   eCM      dsigma/dt_pp        dsigma/dt_ppbar "
       << "\n           ABMST      RPP      ABMST      RPP" << endl;
  for (int iE = 0; iE < 4; ++iE) {
    double eCM        = eCMval[iE];
    double s          = eCM * eCM;
    double sigPPABMST = dsigmadtABMST(s, tNow, true);
    double sigPbABMST = dsigmadtABMST(s, tNow, false);
    double sigPPRPP   = dsigmadtRPP(s, tNow, true);
    double sigPbRPP   = dsigmadtRPP(s, tNow, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(5)
         << setw(10) << sigPPABMST << setw(10) << sigPPRPP
         << setw(10) << sigPbABMST << setw(10) << sigPbRPP << endl;
  }
  */

  /*
  // t dependence of dsigma_el/dt.
  double eCM = 8000.;
  double s   = eCM * eCM;
  cout << "\n    |t|         dsigma/dt_pp            dsigma/dt_ppbar "
       << "\n              ABMST        RPP        ABMST        RPP" << endl;
  for (int iT = 0; iT < nT; ++iT) {
    double t           = tval[iT];
    double dsigPPABMST = dsigmadtABMST( s, t, true);
    double dsigPbABMST = dsigmadtABMST( s, t, false);
    double dsigPPRPP   = dsigmadtRPP( s, t, true);
    double dsigPbRPP   = dsigmadtRPP( s, t, false);
    cout << scientific << setprecision(5) << setw(8) << t
         << setw(12) << dsigPPABMST << setw(12) << dsigPPRPP
         << setw(12) << dsigPbABMST << setw(12) << dsigPbRPP << endl;
  }
  */

  // dsigma_el/dt spectra.
  double eCM = 23.;
  double s   = eCM * eCM;
  Hist dsigHPPABMST( "dsigma/dt pp ABMST", 100, -0.025, 4.975);
  Hist dsigHPbABMST( "dsigma/dt ppbar ABMST", 100, -0.025, 4.975);
  Hist dsigHPPRPP( "dsigma/dt pp RPP", 100, -0.025, 4.975);
  Hist dsigHPbRPP( "dsigma/dt ppbar RPP", 100, -0.025, 4.975);
  for (int iT = 0; iT < 100; ++iT) {
    double t = -0.05 * iT;
    dsigHPPABMST.fill( -t, dsigmadtABMST( s, t, true) );
    dsigHPbABMST.fill( -t, dsigmadtABMST( s, t, false) );
    dsigHPPRPP.fill(   -t, dsigmadtRPP( s, t, true) );
    dsigHPbRPP.fill(   -t, dsigmadtRPP( s, t, false) );
  }
  dsigHPPABMST.takeLog();
  dsigHPbABMST.takeLog();
  dsigHPPRPP.takeLog();
  dsigHPbRPP.takeLog();
  cout << dsigHPPABMST << dsigHPPRPP << dsigHPbABMST << dsigHPbRPP;


  // Comparison with Peter.
  double tval[4] = {-0.0001, -0.001, -0.003, -0.01};
  double sval = 53. * 53.;
  for (int iT = 0; iT < 4; ++iT) cout << scientific << setprecision(4)
    << " t = " << tval[iT] << " dsigma/dt = "
    << dsigmadtABMST( sval, tval[iT], true) << endl;

  return 0;
}
