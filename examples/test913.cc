// test913.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Particle Physics Phenomenology 2015 task 7.
// Emission of gluons off a b quark down to a cutoff scale around m_b.
// Pythia is only used for random numbers and histograms.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Stopping scale of evolution; to be tuned.
  double QMin   = 4.9;

  // Main parameters of the run.
  int    nEvt   = 10000;              // Number of events.
  //double eJet   = 45.;                // Jet energy, part a.
  double eJet   = 1000.;              // Jet energy, part b.
  double mb     = 4.5;                // b quark mass.
  double Lambda = 0.2;                // Lambda in 5-flavour running alpha_s.
  double epsi   = 1e-4;               // Soft-gluon cutoff.

  // Pythia generator, used to extract random number generator.
  Pythia pythia;
  Rndm& rndm = pythia.rndm;

  // Book histograms.
  Hist xbEnd( "b( x, Qmin)",   100, 0., 1.0001);

  // Use t = ln(Q^2 / Lambda^2) as evolution variable.
  double b0   = 23. / 6.;
  double tMax = 2. * log( eJet / Lambda);
  double tMin = 2. * log( QMin / Lambda);
  double zMax = 1. - epsi;
  double xNow, tNow, zNow;

  // Loop over number of evolutions.
  for (int iEvt = 0; iEvt < nEvt; ++ iEvt) {

    // Initially b quark carries the full momentum.
    xNow = 1.;
    tNow = tMax;

    // Evolve towards lower Q scales.
    for ( ; ; ) {

      // Range of possible z values and overestimated branching kernel.
      double zMin = mb / (xNow * eJet);
      if (zMin > zMax) break;
      double intZ = (8./3.) * log((1. - zMin) / (1. - zMax));

      // Evolution in t. Done if reached lower cutoff.
      tNow *= pow( rndm.flat(), b0 / intZ);
      if (tNow < tMin) break;

      // Choice of z and, if accepted, reduce x accordingly.
      zNow = 1. - (1. - zMin) * pow((1. - zMax) / (1. - zMin), rndm.flat());
      if (1. + pow2(zNow) > 2. * rndm.flat()) xNow *= zNow;
    }

    // Fill histogram.
    xbEnd.fill( xNow);
  }

  // Print histogram.
  cout << xbEnd;

  // Done.
  return 0;
}
