// test140.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Hadronization of top using the R-hadron machinery, in e+e- collisions.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// H1 Fit 2 Pomeron flux.
// flux = norm * exp(B_Pom*t)/x^(2*\alpha(t)-1).

double pomFluxH1Fit2(double x, int ifit, bool isReggeon) {

  // Set up constants. Values from fit 2.
  double ap = (isReggeon) ? 0.90 : 0.26;
  double b0 = (isReggeon) ? 2.00 : 4.60;
  double a0 = (isReggeon) ? 0.57 : 1.20;
  double norm = 0.;
  if (ifit == 2) norm = (isReggeon) ? 16. : 1.;
  if (ifit == 3) norm = (isReggeon) ? 15.9 : 0.;
  double flux = 0.;
  if (ifit < 2 || ifit > 3){
    cout << "H1 Fit iFlux=" << ifit << ": Not a valid fit."
         << " Must be either 2 or 3" << endl;
    return 0.;
  }

  // t range from t_max = -mp^2x^2/(1-x) to t_min = -1.
  double mp    = 0.93827231;
  double tMin  = -1.;
  double tMax  = -pow(mp*x, 2.)/(1. - x);
  if (tMax < tMin) return 0.;

  // Flux.
  double b = b0 + 2. * ap * log(1./x);
  flux = norm * exp(log(1/x) * ( 2. * a0 - 1.));
  flux *= (exp(b * tMax) - exp(b * tMin)) / b;

  return flux;
}

//==========================================================================

// Integrated H1 Fit 2 Pomeron flux.
// flux = norm * exp(B_Pom*t)/x^(2*\alpha(t)-1).

double integrateH1Fit2Flux(int iFit, bool isReggeon) {

  // Number of points, x ranges and initial values.
  int    nXLin  = 100;
  double xLow   = 0.035;
  double xUpp   = 0.095;
  double dxLin  = (xUpp - xLow) / nXLin;
  double sum    = 0.;

  // Integral over xPom range.
  for (int iLin = 0; iLin < nXLin; ++iLin) {
    double x = xLow + (iLin + 0.5) * dxLin;

    // Get integral
    double flux = pomFluxH1Fit2(x, iFit, isReggeon);
    sum += dxLin * flux;
  }

  return sum;
}

//==========================================================================

// Integration to check momentum sum rule.
// id = 0 gives sum of all flavours, else specified flavour.

double integrate(PDF* nowPDF, double Q2now, int id = 0) {

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
    if (id == 0) {
      sumNow = nowPDF->xf( 21, x, Q2now);
      for (int i = 1; i < 6; ++i)
        sumNow += nowPDF->xf( i, x, Q2now) + nowPDF->xf( -i, x, Q2now);
    } else sumNow = nowPDF->xf( id, x, Q2now);
    sum   += dxLin * sumNow;
  }

  // Integration at small x in logarithmic steps.
  for (int iLog = 0; iLog < nLog; ++iLog) {
    x      = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
    if (id == 0) {
      sumNow = nowPDF->xf( 21, x, Q2now);
      for (int i = 1; i < 6; ++i)
        sumNow += nowPDF->xf( i, x, Q2now) + nowPDF->xf( -i, x, Q2now);
    } else sumNow = nowPDF->xf( id, x, Q2now);
    sum   += dxLog * x * sumNow;
  }

  // Done.
  return sum;

}

//==========================================================================

int main() {

  // Trial values.
  double Q2 = 75.;
  double x  = 0.1;

  /*
  // Integrated Pomeron flux in 0.035 < xP < 0.095.
  double normFlux = integrateH1Fit2Flux(2, false);
  cout << " normFlux = " << normFlux << endl;

  // Initialize Pomeron PDF.
  PDF* pomB = new PomH1FitAB( 990, 2);

  // Check normalization of Pomeron PDF.
  double normPDF = integrate( pomB, Q2);
  cout << " normPDF  = " << normPDF << endl;

  // Histogram FtildeD spectrum.
  Hist FtildeD( "Ftilde_JJ^D(beta)", 100, 0.00, 1.0);
  for (int iB = 0; iB < 100; ++iB) {
    double x = 0.01 * (iB + 0.5);
    double sum = pomB->xf( 21, x, Q2) + (4./9.) * ( pomB->xf( 1, x, Q2)
      + pomB->xf( 2, x, Q2) + pomB->xf( 3, x, Q2) + pomB->xf( 4, x, Q2)
      + pomB->xf( 5, x, Q2) + pomB->xf( -1, x, Q2) + pomB->xf( -2, x, Q2)
      + pomB->xf( -3, x, Q2) + pomB->xf( -4, x, Q2) + pomB->xf( -5, x, Q2) );
  FtildeD.fill( x, normFlux * sum);
  }

  // Done.
  cout << FtildeD;
  */

  // Check output from ACTW Pomeron PDFs + more for comparison.
  for (int iP = 0; iP < 9; ++iP) {
    cout << "\n Begin PDF no = " << iP << endl;
    PDF* pomACTW;
    if      (iP == 0) pomACTW = new CTEQ6pdf(2212, 1);
    else if (iP < 4)  pomACTW = new PomH1FitAB( 990, iP);
    else if (iP == 4) pomACTW = new PomH1Jets( 990);
    else              pomACTW = new CTEQ6pdf(990, 6 + iP);

    // Check normalization of Pomeron PDF.
    double normPDF = integrate( pomACTW, Q2);
    cout << " normPDF  = " << normPDF << endl;

    // Print values.
    for (int i = -6; i < 7; ++i) {
      int id = (i == 0) ? 21 : i;
      cout << " id = "    << setw(3) << id << fixed << setprecision(5)
           << " xf = "    << setw(9) << pomACTW->xf(id,x,Q2)
           << " xfVal = " << setw(9) <<  pomACTW->xfVal(id,x,Q2)
           << " xfSea = " << setw(9) <<  pomACTW->xfSea(id,x,Q2)
           << " int = "   << setw(9) << integrate( pomACTW, Q2, id) << endl;
    }

    // Done with current PDF set.
    delete pomACTW;
  }

  return 0;
}
