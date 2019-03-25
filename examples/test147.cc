// New parametrizations for sigma_tot and dsigma_el/dt
// by Donnachie-Landshoff and the Review of Particle Physics.

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
const double SPROTON = 0.880;
const double HBARC2  = 0.38938;

// Common method.
complex sModAlp( double sMod, double alpha) {
  return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha);
}

//==========================================================================

// Donnachie-Landshoff 1309.1292. t is now defined negative.

const double EPSIDL[]  = { 0.1103784, -0.3268499, -0.5046476};
const double ALPPDL[]  = { 0.1654554, 0.8, 0.92};
const double NORMDL[]  = { 338.9625, 211.8584, 103.9416};
const double SLOPEDL[] = { 7.854109, 2.470217};
const double FRACSDL[] = { 0.6822898, 0.3177102};
const double TRIGDL[]  = { 3.076837, 4.230092};

//--------------------------------------------------------------------------

// Amplitude.

complex amplitudeDL( double s, double t, double ispp = true) {

  // Common values.
  double snu  = s - 2. * SPROTON + 0.5 * t;
  double ampt = FRACSDL[0] * exp(SLOPEDL[0] * t)
              + FRACSDL[1] * exp(SLOPEDL[1] * t);

  // Single pomeron exchange.
  double alphaPom = 1. + EPSIDL[0] + ALPPDL[0] * t;
  complex ampPom  = -NORMDL[0] * ampt * sModAlp( ALPPDL[0] * snu, alphaPom);

  // Even reggeon exchange.
  double alphaPos = 1. + EPSIDL[1] + ALPPDL[1] * t;
  complex ampPos  = -NORMDL[1] * ampt * sModAlp( ALPPDL[1] * snu, alphaPos);

  // Odd reggeoon exchange.
  double alphaNeg = 1. + EPSIDL[2] + ALPPDL[2] * t;
  complex ampNeg  = complex( 0. ,NORMDL[2]) * ampt
                  * sModAlp( ALPPDL[2] * snu, alphaNeg);

  // Two-pomeron exchange.
  double alphaPP  = 1. + 2. * EPSIDL[0] + 0.5 * ALPPDL[0] * t;
  complex aPlogL  = ALPPDL[0] *  complex( log(ALPPDL[0] * snu), -0.5 * M_PI);
  complex ampPP   = pow2(NORMDL[0]) / (32. * M_PI)
    * sModAlp( ALPPDL[0] * snu, alphaPP) * ALPPDL[0] *
    ( pow2(FRACSDL[0]) * exp( 0.5 * SLOPEDL[0] * t) / (SLOPEDL[0] + aPlogL)
    + pow2(FRACSDL[1]) * exp( 0.5 * SLOPEDL[1] * t) / (SLOPEDL[1] + aPlogL) );

  // Triple-gluon exchange.
  double amp3g    = (t < -TRIGDL[1]) ? TRIGDL[0] / pow4(t)
    :  TRIGDL[0] * exp(2. - 2. * pow2(t / TRIGDL[1])) / pow4(TRIGDL[1]);

  // Add up contributions.
  complex ampSum = (ampPom + ampPos + ((ispp) ? -ampNeg : ampNeg) + ampPP)
    / snu + ((ispp) ? amp3g : -amp3g);
  //cout << fixed << setprecision(5) << " DL  (re, im) = ( "
  //     << HBARC2 * real(ampSum)
  //     << " , " << HBARC2 * imag(ampSum) << " ) " << endl;
  return ampSum;

}

//--------------------------------------------------------------------------

// Total cross section.

double sigmaDL(double s, double ispp = true) {
  return HBARC2 * imag(amplitudeDL( s, 0., ispp));
}

//--------------------------------------------------------------------------

// The rho parameter.

double rhoDL(double s, double ispp = true) {
  complex amp = amplitudeDL( s, 0., ispp);
  return real(amp) / imag(amp);
}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double dsigmadtDL(double s, double t, double ispp = true) {
  return HBARC2 * pow2(abs(amplitudeDL( s, t, ispp))) / (16. * M_PI);
}

//--------------------------------------------------------------------------

// Slope parameter , defined between t = 0 and t = -0.01 by default.

double bSlopeDL(double s, double ispp = true, double t2 = -0.01) {
  return log( dsigmadtDL(s, t2, ispp) / dsigmadtDL(s, 0., ispp) ) / t2;
}

//--------------------------------------------------------------------------

// Total elastic cross section, by integration or optical theorem.

double sigmaelDL(double s, double ispp = true, bool integrate = true) {
  double sigel = 0.;
  if (integrate) {
    int nPoints = 1000;
    double samplingSlope = 10.;
    for (int i = 0; i < nPoints; ++i) {
      double y = (i + 0.5) / nPoints;
      double t = log(y) / samplingSlope;
      sigel += dsigmadtDL(s, t, ispp) / y;
    }
    sigel /= (nPoints * samplingSlope);
  } else {
    sigel = pow2(sigmaDL(s, ispp)) * (1. + pow2(rhoDL(s, ispp)))
      / (HBARC2 * 16. * M_PI * bSlopeDL(s, ispp));
  }
  return sigel;
}

//====================================================================

// Total cross section according to 2014 RPP.

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
  double  shat   = s - 2. * SPROTON - 0.5 * t;
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
  complex fPosRP = t * (complex( 0., NORMRPP[6]) / stlog) * exp(BRPP[6] * t)
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
  //cout << fixed << setprecision(5) << " RPP (re, im) = ( "
  //     << real(ampSum) / sqrt(s * ( s - 4. * SPROTON)) << " , "
  //     << imag(ampSum) / sqrt(s * ( s - 4. * SPROTON)) << " ) " << endl;
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
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double sigPPDL  = sigmaDL(s, true);
    double sigPbDL  = sigmaDL(s, false);
    double sigPPRPP = sigmaRPP(s, true);
    double sigPbRPP = sigmaRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << sigPPDL << setw(10) << sigPPRPP
         << setw(10) << sigPbDL << setw(10) << sigPbRPP << endl;
  }

  // Energy dependence of the rho parameter.
  cout << "\n   eCM         rho_pp             rho_ppbar "
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double rhoPPDL  = rhoDL(s, true);
    double rhoPbDL  = rhoDL(s, false);
    double rhoPPRPP = rhoRPP(s, true);
    double rhoPbRPP = rhoRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << rhoPPDL << setw(10) << rhoPPRPP
         << setw(10) << rhoPbDL << setw(10) << rhoPbRPP << endl;
  }

  // Energy dependence of sigma_elastic; nmerical integration.
  cout << "\n   sigma_elastic by numerical integration"
       << "\n   eCM       sigmael_pp         sigmael_ppbar "
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double sigPPDL  = sigmaelDL(s, true, true);
    double sigPbDL  = sigmaelDL(s, false, true);
    double sigPPRPP = sigmaelRPP(s, true, true);
    double sigPbRPP = sigmaelRPP(s, false, true);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << sigPPDL << setw(10) << sigPPRPP
         << setw(10) << sigPbDL << setw(10) << sigPbRPP << endl;
  }

  // Energy dependence of sigma_elastic; optical theorem + fix slope.
  cout << "\n   sigma_elastic via optical theorem + fix slope"
       << "\n   eCM       sigmael_pp         sigmael_ppbar "
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double sigPPDL  = sigmaelDL(s, true, false);
    double sigPbDL  = sigmaelDL(s, false, false);
    double sigPPRPP = sigmaelRPP(s, true, false);
    double sigPbRPP = sigmaelRPP(s, false, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << sigPPDL << setw(10) << sigPPRPP
         << setw(10) << sigPbDL << setw(10) << sigPbRPP << endl;
  }

  // Approximate t slopes.
  cout << "\n   eCM          B_pp             B_ppbar "
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < nE; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double bPPDL  = bSlopeDL(s, true);
    double bPbDL  = bSlopeDL(s, false);
    double bPPRPP = bSlopeRPP(s, true);
    double bPbRPP = bSlopeRPP(s, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << bPPDL << setw(10) << bPPRPP
         << setw(10) << bPbDL << setw(10) << bPbRPP << endl;
  }

  /*
  // Energy dependence of dsigma/dt at some fixed t; here -0.1.
  double tNow = -0.1;
  cout << "\n   dsigma/dt at t = " << fixed << setprecision(4) << tNow
       << "\n   eCM      dsigma/dt_pp        dsigma/dt_ppbar "
       << "\n            DL       RPP        DL       RPP" << endl;
  for (int iE = 0; iE < 4; ++iE) {
    double eCM      = eCMval[iE];
    double s        = eCM * eCM;
    double sigPPDL  = dsigmadtDL(s, tNow, true);
    double sigPbDL  = dsigmadtDL(s, tNow, false);
    double sigPPRPP = dsigmadtRPP(s, tNow, true);
    double sigPbRPP = dsigmadtRPP(s, tNow, false);
    cout << fixed << setprecision(0) << setw(6) << eCM << setprecision(3)
         << setw(10) << sigPPDL << setw(10) << sigPPRPP
         << setw(10) << sigPbDL << setw(10) << sigPbRPP << endl;
  }
  */

  /*
  // t dependence of dsigma_el/dt.
  double eCM = 8000.;
  double s   = eCM * eCM;
  cout << "\n    |t|         dsigma/dt_pp            dsigma/dt_ppbar "
       << "\n               DL         RPP          DL         RPP" << endl;
  for (int iT = 0; iT < nT; ++iT) {
    double t = tval[iT];
    double dsigPPDL  = dsigmadtDL( s, t, true);
    double dsigPbDL  = dsigmadtDL( s, t, false);
    double dsigPPRPP = dsigmadtRPP( s, t, true);
    double dsigPbRPP = dsigmadtRPP( s, t, false);
    cout << scientific << setprecision(3) << setw(8) << t
         << setw(12) << dsigPPDL << setw(12) << dsigPPRPP
         << setw(12) << dsigPbDL << setw(12) << dsigPbRPP << endl;
  }
  */

  /*
  // dsigma_el/dt spectra.
  double eCM = 23.;
  double s   = eCM * eCM;
  Hist dsigHPPDL( "dsigma/dt pp DL", 100, -0.025, 4.975);
  Hist dsigHPbDL( "dsigma/dt ppbar DL", 100, -0.025, 4.975);
  Hist dsigHPPRPP( "dsigma/dt pp RPP", 100, -0.025, 4.975);
  Hist dsigHPbRPP( "dsigma/dt ppbar RPP", 100, -0.025, 4.975);
  for (int iT = 0; iT < 100; ++iT) {
    double t = -0.05 * iT;
    dsigHPPDL.fill(  -t, dsigmadtDL( s, t, true) );
    dsigHPbDL.fill(  -t, dsigmadtDL( s, t, false) );
    dsigHPPRPP.fill( -t, dsigmadtRPP( s, t, true) );
    dsigHPbRPP.fill( -t, dsigmadtRPP( s, t, false) );
  }
  dsigHPPDL.takeLog();
  dsigHPbDL.takeLog();
  dsigHPPRPP.takeLog();
  dsigHPbRPP.takeLog();
  cout << dsigHPPDL << dsigHPPRPP << dsigHPbDL << dsigHPbRPP;
  */

  // Comparison with Peter.
  double tval[4] = {-0.0001, -0.001, -0.003, -0.01};
  double sval = 53. * 53.;
  for (int iT = 0; iT < 4; ++iT) cout << scientific << setprecision(4)
    << " t = " << tval[iT] << " dsigma/dt = "
    << dsigmadtDL( sval, tval[iT], true) << endl;

  return 0;
}
