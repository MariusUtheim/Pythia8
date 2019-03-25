// test173.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of Bessel K_{1/4} implementation.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Evaluate Bessel function K_{1/4}(x).
// Use power series for x < 2.5 and asymptotic expansion for x > 2.5.
// Number of terms picked to have accuracy better than 1 per mille.
// Based on M. Abramowitz and I.A. Stegun, eqs. 9.6.2, 9.6.10, 9.7.2.

double BesselK14(double x) {

  // Power series expansion of K_{1/4} : k = 0 term.
  if (x < 2.5) {
    double xRat  = 0.25 * x * x;
    double prodP = pow( 0.5 * x, -0.25) / 1.2254167024;
    double prodN = pow( 0.5 * x,  0.25) / 0.9064024771;
    double sum   = prodP - prodN;

    // Power series expansion of K_{1/4} : m > 0 terms.
    for (int k = 1; k < 6; ++k) {
      prodP *= xRat / (k * (k - 0.25));
      prodN *= xRat / (k * (k + 0.25));
      sum   += prodP - prodN;
    }
    sum *= M_PI * sqrt(0.5);
    return sum;

  // Asymptotic expansion of K_{1/4}.
  } else {
    double asym  = sqrt(M_PI * 0.5 / x) * exp(-x);
    double term1 = -         0.75 / ( 8. * x);
    double term2 = -term1 *  8.75 / (16. * x);
    double term3 = -term2 * 24.75 / (24. * x);
    double term4 = -term3 * 48.75 / (32. * x);
    asym *= 1. + term1 + term2 + term3 + term4;
    return asym;
  }
}

//==========================================================================

// Pick an x according to K_{1/4}(x) * x^{3/4}, using upper estimate.

double pickX( double fracSmallX, Rndm& rndm) {
  double x, approx, wanted;
  do {
    x = (rndm.flat() < fracSmallX) ? rndm.flat() : 1. - log(rndm.flat()) / 0.9;
    approx = (x < 1.) ? 0.6 : 1.2 * exp(-0.9 * x);
    wanted = BesselK14(x) * pow( x, 0.75);
  } while (rndm.flat() * approx > wanted);
  return x;
}

//==========================================================================

int main() {

  // Table for values, compared with proposed overestimate.
  double x, approx, wanted, x1, x2;
  for (int ix = 1; ix <= 60; ++ix) {
    x      = (ix < 50) ? 0.1 * ix : ix - 45.;
    approx = (x < 1.) ? 0.6 : 1.2 * exp(-0.9 * x);
    wanted = BesselK14(x) * pow( x, 0.75);
    cout << fixed << setprecision(4) << setw(10) << x << setw(10)
         << wanted << fixed << setw(10) << wanted / approx << endl;
  }

  // Use Pythia to get a random number generator.
  Pythia pythia(" ", false);
  Rndm rndm = pythia.rndm;

  // Histogram correct distribution.
  Hist wadist("wanted K_{1/4}(x) * x^{3/4}", 100, 0., 5.);
  for (int ix = 0; ix < 10000; ++ix) {
    x      = 0.0005 * (ix + 0.5);
    wanted = BesselK14(x) * pow( x, 0.75);
    wadist.fill( x, wanted);
  }

  // Calculate Monte Carlo efficiency for overestimate.
  int nR   = 100000000;
  int nAcc = 0;
  Hist apdist("approximate distribution", 100, 0., 5.);
  double fracSmallX = 0.6 / (0.6 + (1.2/0.9) * exp(-0.9));
  for (int iR = 0; iR < nR; ++iR) {
    x = (rndm.flat() < fracSmallX) ? rndm.flat() : 1. - log(rndm.flat()) / 0.9;
    approx = (x < 1.) ? 0.6 : 1.2 * exp(-0.9 * x);
    wanted = BesselK14(x) * pow( x, 0.75);
    if (wanted > approx) cout << " Error: for x = " << scientific
      << setprecision(4) << x << " ratio is = " << wanted / approx << endl;
    if (rndm.flat() * approx < wanted) {
      ++nAcc;
      apdist.fill( x);
    }
  }
  cout << "\n Efficiency of x selection is " << fixed << setprecision(4)
       << double(nAcc)/double(nR) << endl;

  // Print histograms.
  Hist apwadist("ratio approx/wanted", 100, 0., 5.);
  apwadist = apdist / wadist;
  cout << wadist << apdist << apwadist;

  // Histogram x * exponential.
  Hist xexpdist("wanted x * exp(-x)", 100, 0., 5.);
  for (int ix = 0; ix < 10000; ++ix) {
    x      = 0.0005 * (ix + 0.5);
    xexpdist.fill( x, x * exp(-x));
  }

  // Study convolution to see if it gives correct answer.
  Hist convdist("convoluted distribution", 100, 0., 5.);
  for (int iR = 0; iR < nR; ++iR) {
    x1 = pickX( fracSmallX, rndm);
    x2 = pickX( fracSmallX, rndm);
    x = sqrt( x1 * x1 + x2 * x2 + 2. * x1 * x2 * cos(M_PI * rndm.flat()) );
    convdist.fill( x);
  }

  // Print histograms.
  Hist convxexpdist("ratio convoluted/(x * exp(-x))", 100, 0., 5.);
  convxexpdist = convdist / xexpdist;
  cout << xexpdist << convdist << convxexpdist;

  return 0;
}
