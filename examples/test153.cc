// test153.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF6 interface and whether PDF's behave sensibly.

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

  // Optional parts of tests.
  bool doMomentumSum   = true;
  bool compareInterpol = true;
  bool checkPomeron    = false;

  // Pythia instance to get random numbers.
  Pythia pythia;

  // Pointers to old default and new tryout PDF sets.
  // Compare same PDF with different interpolation.
  PDF* oldPDF = new NNPDF(2212, 2);
  //PDF* newPDF = new LHAGrid1(2212, "2");
  PDF* newPDF = new LHAGrid1(2212, "lhagrid1:2");
  //PDF* newPDF = new LHAGrid1(2212, "NNPDF23_lo_as_0119_qed_0000.dat");
  //PDF* newPDF = new LHAGrid1(2212, "/Users/torbjorn/code/pythia82"
  //  "/pythia8218/share/Pythia8/xmldoc/NNPDF23_lo_as_0119_qed_0000.dat");
  // Compare different pomeron PDFs.
  //PDF* oldPDF = new PomH1FitAB( 990, 0, 1.);
  //PDF* oldPDF = new PomH1Jets( 990, 1, 1.);
  //PDF* newPDF = new LHAGrid1(990, 1);
  //PDF* newPDF = new LHAGrid1(990, 0,
  //  "../share/Pythia8/xmldoc/hf_pdf_0000.dat");

  // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
  // Default behaviour is to freeze PDF's at boundaries.
  newPDF->setExtrapolate(true);

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
  if (doMomentumSum) {
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
  }

  // Comparison of value at random points.
  if (compareInterpol) {
    cout << endl;
    Hist diff("(new - old)/(new + old)", 100, -0.02, 0.02);
    int id[6] = { 21, 1, -1, 2, -2, 3};
    double x, Q2, oldVal, newVal, diffrat;
    for (int iTry = 0; iTry < 1000000; ++iTry) {
      x  = pow( 1e-6, pythia.rndm.flat() );
      Q2 = pow( 1e6,  pythia.rndm.flat() );
      for (int iid = 0 ; iid < 6; ++iid) {
        oldVal  = oldPDF->xf( id[iid], x, Q2);
        newVal  = newPDF->xf( id[iid], x, Q2);
        diffrat = (newVal - oldVal) / (newVal + oldVal);
        if (oldVal > 1e-5 || newVal > 1e-5) {
          diff.fill( diffrat );
          if (abs(diffrat) > 0.02) cout << scientific << " x = " << x
            << " Q2 = " << Q2 << " id = " << id[iid] << " old = " << oldVal
            << " new = " << newVal << " ratio = " <<diffrat << endl;
        }
      }
    }
    cout << diff;

    // Special check at small Q2 and large x.
    Hist uLargexOld( "xu(x) at large x, old", 100, 0.95, 1.00);
    Hist uLargexNew( "xu(x) at large x, new", 100, 0.95, 1.00);
    for (int iX = 0; iX < 100; ++iX) {
      double xLin = 0.95 + 0.05 * 0.01 * (iX + 0.5);
      uLargexOld.fill( xLin, oldPDF->xf( 2, xLin, 1.4) );
      uLargexNew.fill( xLin, newPDF->xf( 2, xLin, 1.4) );
    }
    cout << uLargexOld << uLargexNew << uLargexNew / uLargexOld;
  }

  // Shape of new PDF in linear x scale.
  if (checkPomeron) {
    Hist newg("new g PDF, lin x", 100, 0., 1.);
    Hist newq("new light q PDF, lin x", 100, 0., 1.);
    Hist newc("new c PDF, lin x", 100, 0., 1.);
    Hist newb("new b PDF, lin x", 100, 0., 1.);
    double Q2Lin = 100.;
    double xLin, updf, qpdf;
    int idchk[5] = { 1, 3, -1, -2, -3};
    for (int iX = 0; iX < 100; ++iX) {
      xLin = 0.01 * (iX + 0.5);
      newg.fill( xLin, newPDF->xf( 21, xLin, Q2Lin) );
      updf = newPDF->xf( 2, xLin, Q2Lin);
      newq.fill ( xLin, updf);
      for (int i = 0; i < 5; ++i) {
        qpdf = newPDF->xf( idchk[i], xLin, Q2Lin);
        if ( abs(qpdf - updf) > 1e-6 * updf)
          cout << " Error: mismatch at x = " << xLin << " for id = "
               << idchk[i] << " of pdf = " << scientific << setprecision(3)
               << qpdf << " vs expected = " << updf << endl;
      }
      newc.fill( xLin, newPDF->xf( 4, xLin, Q2Lin) );
      newb.fill( xLin, newPDF->xf( 5, xLin, Q2Lin) );
    }
    cout << newg << newq << newc << newb;
  }

  // Done.
  delete oldPDF;
  delete newPDF;
  return 0;
}
