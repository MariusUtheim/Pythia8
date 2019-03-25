// test170.cc is a part of the PYTHIA event generator.
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
  pythia.readString("Higgs:useBSM = on");
  pythia.readString("HiggsBSM:gg2H2ttbar = on");
  pythia.readString("HiggsBSM:qqbar2H2ttbar = on");

  // Force decay chain H -> A A -> 4 gamma.
  pythia.readString("35:m0 = 750.");
  pythia.readString("35:mWidth = 0.00407");
  pythia.readString("35:doForceWidth = on");
  pythia.readString("35:onMode = off");
  pythia.readString("35:onIfMatch = 36 36");
  pythia.readString("36:m0 = 1.1");
  pythia.readString("36:mMin = 1.0");
  pythia.readString("36:mWidth = 0.001");
  pythia.readString("36:doForceWidth = on");
  pythia.readString("36:onMode = off");
  pythia.readString("36:onIfMatch = 22 22");

  // Initialize.
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
