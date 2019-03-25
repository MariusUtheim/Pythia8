// test918.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Particle Physics Phenomenology reexam 2018 task 7.
// Toy PDF evolution for Particle Physics Phenomenology course exam.
// Pythia is only used for random numbers and histograms.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Main parameters of the run.
  int    nEvt   = 100000;              // Number of events.
  double Q0     = 0.5;                // Starting scale of evolution.
  double QMax   = 10.;                // Finishing scale of evolution.
  double alphaS = 0.15;               // Fixed alpha_strong value.
  double epsi   = 1e-4;               // Soft-gluon cutoff.
  double xNow, zNow, QNow;

  // Pythia generator, used to extract random number generator.
  Pythia pythia;
  Rndm& rndm = pythia.rndm;

  // Book histograms.
  Hist uBeg( "u( x, 0.5 GeV)", 100, 0., 1.);
  Hist uEnd( "u( x, 10 GeV)" , 100, 0., 1.);
  Hist xuBeg( "xu( x, 0.5 GeV)", 100, 0., 1.);
  Hist xuEnd( "xu( x, 10 GeV)" , 100, 0., 1.);

  // Integral over z for evolution (for overestimation).
  double alpS2pi = alphaS / (2. * M_PI);
  double zMin = 0.;
  double zMax = 1. - epsi;
  double zInt = alpS2pi * (4./3.) * 2. * log((1. - zMin)/(1. - zMax));
  cout << "zInt = " << fixed << setprecision(3) << zInt << endl;

  // Loop over number of evolutions.
  for (int iEvt = 0; iEvt < nEvt; ++ iEvt) {

    // Pick initial x values of quark. Fill histograms.
    QNow = Q0;
    do xNow = rndm.flat();
    while (pow2(xNow) + pow2(1. - xNow) < rndm.flat());
    uBeg.fill( xNow);
    xuBeg.fill( xNow, xNow);

    // Loop over branchings.
    do {

      // Take a step up in Q. Finished if passed QMax.
      do {
        QNow /= pow( rndm.flat(), 0.5 / zInt);
        if (QNow > QMax) break;

        // Pick z and correct for overestimated branching kernel.
        zNow = 1. - (1. - zMin) * pow( (1. - zMax) / (1. - zMin), rndm.flat() );
      } while ( 1. + zNow*zNow < 2. * rndm.flat());

      // Update x value and continue looping until passed QMax.
      if (QNow < QMax) xNow *= zNow;
    } while (QNow < QMax);


    // Fill histograms.
    uEnd.fill( xNow);
    xuEnd.fill( xNow, xNow);

  // End of loop over events.
  }

  // Print histograms.
  cout << uBeg << uEnd << xuBeg << xuEnd;

  // Done.
  return 0;
}
