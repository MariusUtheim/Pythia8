// test230.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Compare partonic and hadronic multiplicities

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 10000;
  int nAbort  = 5;

  // Book histograms.
  Hist nqgH( "hadronizing partonic multiplicity", 100, -0.5, 499.5);
  Hist n71H( "hadronizing partonic multiplicity - 2", 100, -0.5, 499.5);
  Hist nprH( "primary hadronic multiplicity", 100, -0.5, 499.5);

  // Generator setup.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("SoftQCD:nondiffractive = on");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("HadronLevel:Decay = off");
  pythia.init();

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Count number of partons before hadronization.
    int nqg = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) ++nqg;

    // Hadronize, but do not decay, and count primary hadrons.
    pythia.forceHadronLevel();
    if (iEvent == 0) event.list();
    int npr = 0;
    int n71 = 0;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].status() == -71) ++n71;
      if (event[i].isFinal()) ++npr;
    }

    // Fill histograms.
    nqgH.fill( nqg);
    n71H.fill( n71);
    nprH.fill( npr);

  // End of event loop. Print histograms.
  }
  pythia.stat();
  cout << nqgH << n71H << nprH;

  // Done.
  return 0;
}
