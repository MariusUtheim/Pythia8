// test180.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Overlap factors.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Sampling.
  double deltaB = 0.01;
  double normPi = 1. / (2. * M_PI);
  int    nb     = 1000;
  double k      = 33.745;

  // Integrals.
  double Oint  = 0.;
  double Pint  = 0.;
  double OOint = 0.;
  double OPint = 0.;
  double PPint = 0.;

  // Step through b space.
  double bNow = -0.5 * deltaB;
  double bArea = 0.;
  double Onow, Pnow;
  for (int ib = 0; ib < nb; ++ib) {
    bNow += deltaB;
    bArea = 2. * M_PI * bNow * deltaB;

    // Overlap for single Gaussian.
    Onow = normPi * exp(- bNow * bNow);
    Pnow = 1. - exp(- M_PI * k * Onow);

    // Update integrals.
    Oint  += bArea * Onow;
    Pint  += bArea * Pnow;
    OOint += bArea * Onow * Onow;
    OPint += bArea * Onow * Pnow;
    PPint += bArea * Pnow * Pnow;
  }

  // Results.
  cout << fixed << setprecision(3) << " O = " << Oint << " P = " << Pint
       << " OO = " << OOint << " OP = " << OPint << " PP = " << PPint << endl;
  double ans1 = (OOint * Pint) / (OPint * Oint);
  double ans2 = (OOint / PPint) * pow2(Pint / Oint);
  double ans3 = OOint / pow2(Oint);
  double ans4 = (OOint * Pint) / pow2(Oint);
  cout << " answer is " << ans1 << " or " << ans2  << " or " << ans3
       << " or " << ans4 << " ?" << endl;
  double nAvg = M_PI * k * Oint / Pint;
  cout << " nAvg = " << nAvg << endl;

  // Done.
  return 0;

}

// Printout in MultipleInteractions.cc after line 2170:

//   cout << fixed << setprecision(3) << " kNow = " << kNow
//        << " Oint = " << overlapInt << " Pint = " << probInt
//        << " OPint = " << probOverlapInt << endl;
//   cout << " enhancement factor = " << enhanceBinit << endl;
