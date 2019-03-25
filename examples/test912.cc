// test912.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Particle Physics Phenomenology 2015 task 8.
// Toy PDF evolution for Particle Physics Phenomenology course exam.
// Pythia is only used for random numbers and histograms.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Main parameters of the run.
  int    nEvt   = 10000;              // Number of events.
  double Q0     = 1.;                 // Starting scale of evolution.
  double QMax   = 100.;               // Finishing scale of evolution.
  double alphaS = 0.2;                // Fixed alpha_strong value.
  double xMin   = 0.01;               // Do not trace low-x gluons.
  double epsi   = 1e-4;               // Soft-gluon cutoff.

  // Pythia generator, used to extract random number generator.
  Pythia pythia;
  Rndm& rndm = pythia.rndm;

  // Book histograms.
  Hist xgBeg( "xg( x, 1 GeV)",   100, 0., 1.);
  Hist xgEnd( "xg( x, 100 GeV)", 100, 0., 1.);
  Hist ngEnd( "number of gluons with x > 0.01", 100, -0.5, 99.5);

  // Integral over z for evolution (for overestimation).
  double alpS2pi = alphaS / (2. * M_PI);
  double zMin = epsi;
  double zMax = 1. - epsi;
  double zInt = 6 * alpS2pi * log((1. - zMin)/(1. - zMax));

  // Loop over number of evolutions.
  for (int iEvt = 0; iEvt < nEvt; ++ iEvt) {

    // Create vectors of x and Q values.
    vector<double> xVal, QVal;

    // Pick initial x values of two gluons.
    double xBeg = rndm.flat();
    xVal.push_back(xBeg);
    xVal.push_back(1. - xBeg);
    QVal.push_back(Q0);
    QVal.push_back(Q0);

    // Go through list of partons that (maybe) should be evolved.
    int i = -1;
    double QNow, zNow;
    while ( ++i < int(xVal.size()) ) {

      // Do not evolve partons below x = 0.01.
      if (xVal[i] < xMin) continue;

      // Take a step up in Q. Finished if passed QMax.
      QNow = QVal[i];
      do {
        QNow /= pow( rndm.flat(), 0.5/zInt);
        if (QNow > QMax) break;

        // Pick z and correct for overestimated branching kernel.
        zNow = 1. - zMax * pow( zMin / zMax, rndm.flat() );
      } while ( pow2(1. - zNow * (1. - zNow)) < rndm.flat());

      // Bookkeep newly created partons, if any.
      if ( QNow < QMax) {
        xVal.push_back(zNow * xVal[i]);
        xVal.push_back((1. - zNow) * xVal[i]);
        QVal.push_back(QNow);
        QVal.push_back(QNow);

        // Mark branched parton by negating its Q scale.
        QVal[i] = -QVal[i];
      }

    // End loop over partons that may branch.
    }

    // Fill histograms.
    xgBeg.fill( xVal[0], xVal[0]);
    xgBeg.fill( xVal[1], xVal[1]);
    int ngFin = 0;
    for ( i = 0; i < int(xVal.size()); ++i) if (QVal[i] > 0.) {
      xgEnd.fill( xVal[i], xVal[i]);
      if (xVal[i] > xMin) ++ngFin;
    }
    ngEnd.fill( ngFin);

  // End of loop over events.
  }

  // Print histograms.
  cout << xgBeg << xgEnd << ngEnd;

  // Done.
  return 0;
}
