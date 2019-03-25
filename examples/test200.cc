// New parametrizations for sigma_tot and dsigma_el/dt
// by Appleby, Barlow, Molson, Serluca and Toader (ABMST),
// and by the Review of Particle Physics (RPP).

//==========================================================================

#include "Pythia8/Basics.h"
#include "Pythia8/PythiaComplex.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Common values.
const int    NPOINTS    = 1000;
const double SPROTON    = 0.8803544;
const double HBARC2     = 0.38938;
const double ALPHAEM    = 0.00729353;
const double GAMMAEUL = 0.577215665;
const double MINSLOPE   = 10.;

//==========================================================================

// The SigmaABMST class.

// Total and elastic cross section according to
// Appleby, Barlow, Molson, Serluca and Toader (ABMST).

//--------------------------------------------------------------------------

const int    SigmaABMST::NPOINTS  = 1000;
const double SigmaABMST::SPROTON  = 0.8803544;
const double SigmaABMST::HBARC2   = 0.38938;
const double SigmaABMST::ALPHAEM  = 0.00729353;
const double SigmaABMST::GAMMAEUL = 0.577215665;
const double SigmaABMST::MINSLOPE = 10.;
const double SigmaABMST::EPSI[]  = { 0.106231, 0.0972043, -0.510662, -0.302082};
const double SigmaABMST::ALPP[]  = { 0.0449029, 0.278037, 0.821595, 0.904556};
const double SigmaABMST::NORM[]  = { 228.359, 193.811, 518.686, 10.7843};
const double SigmaABMST::SLOPE[] = { 8.38, 3.78, 1.36};
const double SigmaABMST::FRACS[] = { 0.26, 0.56, 0.18};
const double SigmaABMST::TRIG[]  = { 0.3, 5.03};
const double SigmaABMST::LAM2P   = 0.521223;
const double SigmaABMST::BAPPR[] = { 8.5, 1.086};
const double SigmaABMST::LAM2FF  = 0.71;

//--------------------------------------------------------------------------

// Total cross section.

double SigmaABMST::sigmaTot(double s, double ispp = true) {
  return HBARC2 * imag(amplitudeABMST( s, 0., ispp, false));
}

//--------------------------------------------------------------------------

// The rho parameter.

double SigmaABMST::rho(double s, double ispp = true) {
  complex amp = amplitudeABMST( s, 0., ispp, false);
  return real(amp) / imag(amp);
}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double SigmaABMST::dsigmaELdt(double s, double t, double ispp = true,
  bool useCoulomb = false) {
  return HBARC2 * pow2(abs(amplitudeABMST( s, t, ispp, useCoulomb)))
         / (16. * M_PI);
}

//--------------------------------------------------------------------------

// Total elastic cross section, by integration.

double SigmaABMST::sigmaEL(double s, double ispp = true) {

  double sigel = 0.;
  int nPoints = 1000;
  for (int i = 0; i < nPoints; ++i) {
    double y = (i + 0.5) / nPoints;
    double t = log(y) / MINSLOPE;
    sigel += dsigmadtABMST(s, t, ispp, false) / y;
  }
  sigel /= (nPoints * MINSLOPE);
  return sigel;

}

//--------------------------------------------------------------------------

// Select a t value for elastic scattering.

double SigmaABMST::pickELt(double s, Rndm& rndm, double ispp = true,
  bool useCoulomb = false, double tAbsMin = 0.) {

  // Set up sampling details.
  double tRef    = -max( 0., tAbsMin);
  double slope1  = 10.;
  double slope2  = 1.;
  double frac2   = 0.1;
  double rel2    = exp((slope2 - slope1) * tRef) * frac2 / (1. - frac2);
  double sigRef1 = dsigmadtABMST( s, tRef, ispp, false);
  double sigRef2 = dsigmadtABMST( s, tRef - 0.2, ispp, false);
  double sigRef  = (sigRef1 > 2. * sigRef2) ? 2. * sigRef1 : 5. * sigRef2;
  double norm1   = sigRef / (slope1 + rel2 * slope2);
  double norm2   = norm1 * rel2;
  double norm3   = useCoulomb ? -2. * HBARC2 * 4. * M_PI * pow2(ALPHAEM)
                 / tRef : 0.;
  double normTot = norm1 + norm2 + norm3;

  // Repeated tries until accepted.
  double rNow, bNow, tNow, sigNow, sigEst;
  int nTry = 0;
  do {
    ++ nTry;
    rNow = rndm.flat() * normTot;
    if (useCoulomb && rNow > norm1 + norm2) tNow = tRef / rndm.flat();
    else {
      bNow = (rNow < norm1) ? slope1 : slope2;
      tNow = log( rndm.flat() ) / bNow + tRef;
    }
    sigNow = dsigmadtABMST( s, tNow, ispp, useCoulomb);
    sigEst = norm1 * slope1 * exp( slope1 * (tNow-tRef))
           + norm2 * slope2 * exp( slope2 * (tNow-tRef));
    if (useCoulomb) sigEst += norm3 * (-tRef) / pow2(tNow);
    if (sigNow > sigEst) cout << " error ABMST: ratio = " << sigNow/sigEst
      << " for t = " << tNow << (ispp ? " in pp" : " in ppbar") << endl;
  } while (sigNow < sigEst * rndm.flat());
  nTries.fill( nTry );

  return tNow;
}

//--------------------------------------------------------------------------

// Amplitude.

complex SigmaABMST::amplitude( double s, double t, double ispp = true,
  bool useCoulomb = false) {

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
  amp[4]        *= LAM2P * complex( 0., 1.) / (16. * M_PI * snu);

  // Triple-gluon exchange.
  amp[5] = sqrt(16. * M_PI / HBARC2) * TRIG[0] * ((t < -TRIG[1])
         ? 1. / pow4(t) :  exp(4. + 4. * t / TRIG[1]) / pow4(TRIG[1]));

  // Add up contributions.
  complex ampSum = (amp[0] + amp[1] + amp[2] + ((ispp) ? -amp[3] : amp[3])
    + amp[4]) / snu + ((ispp) ? amp[5] : -amp[5]);

  // Optional Coulomb term. Must not be used for t = 0.
  complex ampCou = 0.;
  if (useCoulomb && t < 0.) {
    double bApp  = BAPPR[0] + BAPPR[1] * 0.5 * log(s);
    double phase = (GAMMAEUL  + log( -0.5 * t * (bApp + 8. / LAM2FF)) +
                 - 4. * t / LAM2FF * log(- 4. * t / LAM2FF)
                 - 2. * t / LAM2FF) * (ispp ? 1. : -1.);
    ampCou       = exp( complex( 0., ALPHAEM * phase) ) * 8. * M_PI
                 * ALPHAEM * ampt / t;
  }

  // Combine and return.
  return ispp ? ampSum - ampCou : ampSum + ampCou;

}

//--------------------------------------------------------------------------

// Common method.

complex SigmaABMST::sModAlp( double sMod, double alpha) {
  return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha);
}

//====================================================================

// The SigmaRPP class.

// Total and elastic cross section according to 2014 RPP.

//--------------------------------------------------------------------------

// Common values.
const int    NPOINTS     = 1000;
const double SPROTON     = 0.8803544;
const double HBARC2      = 0.38938;
const double ALPHAEM     = 0.00729353;
const double GAMMAEUL    = 0.577215665;
const double MINSLOPE    = 10.;
const double EPS1RPP[]  = { 1., 0.614, 0.444, 1., 1., 1.};
const double ALPPRPP[]  = { 0.151, 0.8, 0.8, 0.947};
const double NORMRPP[] = { 0.2478, 0.0078, 11.22, -0.150, 148.4, -26.6,
  -1.5, -0.0441, 0., 0.686, -3.82, -8.60, 64.1, 99.1, -58.0, 9.5};
const double BRPP[]    = { 3.592, 0.622, 5.44, 0.205, 5.643, 1.92, 0.41,
  0., 0., 3.013, 2.572, 12.25, 2.611, 11.28, 1.27};
const double KRPP[]     = { 0.3076, 0.0998, 1.678, 0.190, -26.1};

//--------------------------------------------------------------------------

// Total cross section, excluding Coulomb part.

double sigmaRPP(double s, double ispp = true) {
  return imag(amplitudeRPP( s, 0., ispp, false))
         / sqrt(s * ( s - 4. * SPROTON));
}

//--------------------------------------------------------------------------

// The rho parameter.

double rhoRPP(double s, double ispp = true) {
  complex amp = amplitudeRPP( s, 0., ispp, false);
  return real(amp) / imag(amp);
}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double dsigmadtRPP(double s, double t, double ispp = true,
  bool useCoulomb = false) {
  return pow2(abs(amplitudeRPP( s, t, ispp, useCoulomb)))
         / (16. * M_PI * HBARC2 * s * ( s - 4. * SPROTON));
}

//--------------------------------------------------------------------------

// Total elastic cross section, by integration or optical theorem.

double sigmaelRPP(double s, double ispp = true, bool integrate = true) {
  double sigel = 0.;
  if (integrate) {
    int nPoints = 1000;
    for (int i = 0; i < nPoints; ++i) {
      double y = (i + 0.5) / nPoints;
      double t = log(y) / MINSLOPE;
      sigel += dsigmadtRPP(s, t, ispp, false) / y;
    }
    sigel /= (nPoints * MINSLOPE);
  } else {
    sigel = pow2(sigmaRPP(s, ispp)) * (1. + pow2(rhoRPP(s, ispp)))
          / (16. * M_PI * HBARC2 * bSlopeRPP(s, ispp));
  }
  return sigel;
}

//--------------------------------------------------------------------------

// Select a t value for elastic scattering.

double picktRPP(double s, Rndm& rndm, double ispp = true,
  bool useCoulomb = false, double tAbsMin = 0.) {

  // Set up sampling details.
  double tRef    = -max( 0., tAbsMin);
  double slope1  = 10.;
  double slope2  = 1.;
  double frac2   = 0.1;
  double rel2    = exp((slope2 - slope1) * tRef) * frac2 / (1. - frac2);
  double sigRef1 = dsigmadtRPP( s, tRef, ispp, false);
  double sigRef2 = dsigmadtRPP( s, tRef - 0.2, ispp, false);
  double sigRef  = (sigRef1 > 2. * sigRef2) ? 2. * sigRef1 : 5. * sigRef2;
  double norm1   = sigRef / (slope1 + rel2 * slope2);
  double norm2   = norm1 * rel2;
  double norm3   = useCoulomb ? -2. * HBARC2 * 4. * M_PI * pow2(ALPHAEM)
                 / tRef : 0.;
  double normTot = norm1 + norm2 + norm3;

  // Repeated tries until accepted.
  double rNow, bNow, tNow, sigNow, sigEst;
  int nTry = 0;
  do {
    ++ nTry;
    rNow = rndm.flat() * normTot;
    if (useCoulomb && rNow > norm1 + norm2) tNow = tRef / rndm.flat();
    else {
      bNow = (rNow < norm1) ? slope1 : slope2;
      tNow = log( rndm.flat() ) / bNow + tRef;
    }
    sigNow = dsigmadtRPP( s, tNow, ispp, useCoulomb);
    sigEst = norm1 * slope1 * exp( slope1 * (tNow-tRef))
           + norm2 * slope2 * exp( slope2 * (tNow-tRef));
    if (useCoulomb) sigEst += norm3 * (-tRef) / pow2(tNow);
    if (sigNow > sigEst) cout << " error  RPP : ratio = " << sigNow/sigEst
      << " for t = " << tNow << (ispp ? " in pp" : " in ppbar") << endl;
  } while (sigNow < sigEst * rndm.flat());
  nTries.fill( nTry );

  return tNow;
}

//--------------------------------------------------------------------------

// Amplitude.

complex amplitudeRPP(double s, double t, double ispp = true,
  bool useCoulomb = false) {

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

  // Combine nuclear part.
  complex ampSum = ispp ? fPos + fNeg : fPos - fNeg;

  // Optional Coulomb term. Must not be used for t = 0.
  complex ampCou = 0.;
  if (useCoulomb && t < 0.) {
    double bAppr = imag(ampSum) / ( sqrt(s * ( s - 4. * SPROTON))
      * 4. * M_PI * HBARC2 );
    double phase = (log( -0.5 * t * (bAppr + 8. / LAM2FF)) + GAMMAEUL
                 - 4. * t / LAM2FF * log(- 4. * t / LAM2FF)
                 - 2. * t / LAM2FF) * (ispp ? 1. : -1.);
    ampCou       = exp( complex( 0., ALPHAEM * phase) ) * 8. * M_PI * HBARC2
                 * ALPHAEM * s / t * pow(1 - t / LAM2FF, -4.);
  }

  // Combine and return.
  return ispp ? ampSum - ampCou : ampSum + ampCou;

}

//--------------------------------------------------------------------------

// Common method.

complex sModAlp( double sMod, double alpha) {
  return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha);
}

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

//==========================================================================
