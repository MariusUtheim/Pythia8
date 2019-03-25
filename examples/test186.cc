// thest186.cc (main51.cc) is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF interface and whether PDF's behave sensibly.

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

  // Pointers to old default and new tryout PDF sets.
  // LO with high and low alphas, NLO, NNLO in this order.
  Info info;
  PDF* oldPDF = new LHAGrid1(2212, "21");
  PDF* newPDF = new LHAGrid1(2212, "22");
  //PDF* oldPDF = new NNPDF(2212, 1);
  //PDF* newPDF = new LHAGrid1(2212, "17");
  //PDF* oldPDF = new NNPDF(2212, 2);
  //PDF* newPDF = new LHAGrid1(2212, "18");
  //PDF* oldPDF = new NNPDF(2212, 3);
  //PDF* newPDF = new LHAGrid1(2212, "19");
  //PDF* oldPDF = new NNPDF(2212, 4);
  //PDF* newPDF = new LHAGrid1(2212, "20");

  // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
  // Default behaviour is to freeze PDF's at boundaries.
  //newPDF->setExtrapolate(true);

  // Histogram F(x, Q2) = (9/4) x*g(x, Q2) + sum_{i = q, qbar} x*f_i(x, Q2)
  // for range 10^{-8} < x < 1 logarithmic in x and for Q2 = 4 and 100.
  // Also histogram q_val, g and q_sea separately at Q2 = 4..
  Hist oldF4("F( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newF4("F( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratF4("F( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldQV4("x*q_val( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newQV4("x*q_val( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratQV4("q_val( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldG4( "x*g( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newG4( "x*g( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratG4("g( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldQS4("q_sea( x, Q2 = 4) new", 80 , -8., 0.);
  Hist newQS4("q_sea( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratQS4("q_sea( x, Q2 = 4) new/old", 80 , -8., 0.);
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

      // Evaluate old partial PDFs.
      if (iQ == 0) {
        if (x > 1e-5) oldQV4.fill( xLog,
            oldPDF->xf( 1, x, Q2) + oldPDF->xf( 2, x, Q2)
          - oldPDF->xf( -1, x, Q2) - oldPDF->xf( -2, x, Q2) );
        oldG4.fill( xLog, oldPDF->xf( 21, x, Q2) );
        oldQS4.fill( xLog, 2. * ( oldPDF->xf( -1, x, Q2)
          + oldPDF->xf( -2, x, Q2) + oldPDF->xf( -3, x, Q2)
          + oldPDF->xf( -4, x, Q2) + oldPDF->xf( -5, x, Q2) ) );
      }

      // Evaluate new summed PDF.
      double newSum = 2.25 * newPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        newSum += newPDF->xf( i, x, Q2) + newPDF->xf( -i, x, Q2);
      if (iQ == 0) newF4.fill ( xLog, newSum );
      else       newF100.fill ( xLog, newSum );

      // Evaluate new partial PDFs.
      if (iQ == 0) {
         if (x > 1e-5) newQV4.fill( xLog,
            newPDF->xf( 1, x, Q2) + newPDF->xf( 2, x, Q2)
          - newPDF->xf( -1, x, Q2) - newPDF->xf( -2, x, Q2) );
        newG4.fill( xLog, newPDF->xf( 21, x, Q2) );
        newQS4.fill( xLog, 2. * ( newPDF->xf( -1, x, Q2)
          + newPDF->xf( -2, x, Q2) + newPDF->xf( -3, x, Q2)
          + newPDF->xf( -4, x, Q2) + newPDF->xf( -5, x, Q2) ) );
      }

    // End loops over x and Q2 values.
    }
  }

  // Show F(x, Q2), q_val, g and q_sea and their ratio new/old.
  ratF4 = newF4 / oldF4;
  ratQV4 = newQV4 / oldQV4;
  ratG4 = newG4 / oldG4;
  ratQS4 = newQS4 / oldQS4;
  ratF100 = newF100 / oldF100;
  cout << oldF4 << newF4 << ratF4 << oldQV4 << newQV4 << ratQV4
       << oldG4 << newG4 << ratG4 << oldQS4 << newQS4 << ratQS4
       << oldF100 << newF100 << ratF100;

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
