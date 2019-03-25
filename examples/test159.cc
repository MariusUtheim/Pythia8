// test159.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// gamma + gamma -> H -> gamma + gamma.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Energy and process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HiggsSM:gmgm2H = on");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 22 22");

  // Force heavy Higgs to be narrow.
  pythia.readString("25:m0 = 750.");
  pythia.readString("25:mWidth = 1.");
  pythia.readString("25:doForceWidth = on");

  // No event printout.
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize.
  pythia.init();

  // Set up anti-kT jet finder with R = 0.5, pT > 30., |eta| < 5.
  SlowJet slowJet( -1, 0.5, 30., 5.);

  // Histogram.
  Hist massH( "mass of Higgs", 100, 0., 1000.);
  Hist pTj1( "pT jet 1", 100, 0., 200.);
  Hist pTj2( "pT jet 2", 100, 0., 200.);
  Hist yj1( "y jet 1", 100, -5., 5.);
  Hist yj2( "y jet 2", 100, -5., 5.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze event. Remove gammas from H before jet finding.
    massH.fill( pythia.info.mHat() );
    for (int i = 0; i < pythia.event.size(); ++i)
       if (pythia.event[i].id() == 22
       && pythia.event[ pythia.event[i].mother1() ].id() == 25)
       pythia.event[i].statusNeg();
    slowJet.analyze( pythia.event );
    if (slowJet.sizeJet() > 0) {
      pTj1.fill( slowJet.pT(0) );
      yj1.fill( slowJet.y(0) );
    }
    if (slowJet.sizeJet() > 1) {
      pTj2.fill( slowJet.pT(1) );
      yj2.fill( slowJet.y(1) );
    }

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << massH << pTj1 << pTj2 << yj1 << yj2;

  // Done.
  return 0;
}
