// test151.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test extrapolation below xMin.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Integration to check momentum sum rule.

double integrate(PDF* nowPDF, double Q2) {

  // Number of points, x ranges and initial values.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double sum   = 0.;
  double x, sumNow;

  // Integration at large x in linear steps.
  for (int iLin = 0; iLin < nLin; ++iLin) {
    x      = xLin + (iLin + 0.5) * dxLin;
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLin * sumNow;
  }

  // Integration at small x in logarithmic steps.
  for (int iLog = 0; iLog < nLog; ++iLog) {
    x      = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLog * x * sumNow;
  }

  // Done.
  return sum;

}

//==========================================================================

int main() {

  // The Pythia class itself is not used, but some facilities that come along.
  //Pythia pythia;

  // Study Pomeron PDF, frozen at small x or extrapolated.
  //PDF* oldPDF = new CTEQ6pdf( 990, 11);
  //PDF* newPDF = new CTEQ6pdf( 990, 11);
  PDF* oldPDF = new PomH1FitAB( 990, 0);
  PDF* newPDF = new PomH1FitAB( 990, 0);
  //PDF* oldPDF = new PomH1Jets( 990);
  //PDF* newPDF = new PomH1Jets( 990);
  newPDF->setExtrapolate(true);

  // Evaluate one point in phase space.
  cout << scientific << setprecision(3);
  double xA  = 1e-6;
  double Q2A = 100.;
  double pOldA = oldPDF->xf( 21, xA, Q2A);
  double pNewA = newPDF->xf( 21, xA, Q2A);
  cout << " old = " << pOldA << " new = " << pNewA << endl;

  // Histogram F(x, Q2) = (9/4) x*g(x, Q2) + sum_{i = q, qbar} x*f_i(x, Q2)
  // for range 10^{-8} < x < 1 logarithmic in x and for Q2 = 4 and 100.
  Hist oldF4("F( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newF4("F( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratF4("F( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldF100("F( x, Q2 = 100) old", 80 , -8., 0.);
  Hist newF100("F( x, Q2 = 100) new", 80 , -8., 0.);
  Hist ratF100("F( x, Q2 = 100) new/old", 80 , -8., 0.);

  // Loop over the two Q2 values.
  for (int iQ = 0; iQ < 2; ++iQ) {
    double Q2 = (iQ == 0) ? 4. : 100;

    // Loop over x values, in a logarithmic scale
    for (int iX = 0; iX < 80; ++iX) {
      double xLog = -(0.1 * iX + 0.05);
      double x = pow( 10., xLog);

      // Evaluate old summed PDF.
      double oldSum = 2.25 * oldPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        oldSum += oldPDF->xf( i, x, Q2) + oldPDF->xf( -i, x, Q2);
      if (iQ == 0) oldF4.fill ( xLog, oldSum );
      else       oldF100.fill ( xLog, oldSum );

      // Evaluate new summed PDF.
      double newSum = 2.25 * newPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        newSum += newPDF->xf( i, x, Q2) + newPDF->xf( -i, x, Q2);
      if (iQ == 0) newF4.fill ( xLog, newSum );
      else       newF100.fill ( xLog, newSum );

    //End loops over x and Q2 values.
    }
  }

  // Show F(x, Q2) and their ratio new/old.
  ratF4 = newF4 / oldF4;
  ratF100 = newF100 / oldF100;
  cout << oldF4 << newF4 << ratF4 << oldF100 << newF100 << ratF100;

  // Histogram momentum sum as a function of Q2 (or rather log10(Q2)).
  Hist oldXSum("momentum sum(log10(Q2)) old", 100, -2., 8.);
  Hist newXSum("momentum sum(log10(Q2)) new", 100, -2., 8.);

  // Loop over Q2 values.
  for (int iQ = 0; iQ < 100; ++iQ) {
    double log10Q2 = -2.0 + 0.1 * iQ + 0.05;
    double Q2 = pow( 10., log10Q2);

    // Evaluate old and new momentum sums.
    double oldSum = integrate( oldPDF, Q2);
    oldXSum.fill( log10Q2, oldSum);
    double newSum = integrate( newPDF, Q2);
    newXSum.fill( log10Q2, newSum);
  }

  // Show momentum sum as a function of Q2.
  cout << oldXSum << newXSum;

  // Done.
  delete oldPDF;
  delete newPDF;
  return 0;
}
