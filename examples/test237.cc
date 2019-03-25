// test237.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program, illustrating gamma*/Z^0 -> tau'+ tau'-.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Process selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
  pythia.readString("Zprime:gmZmode = 4");
  pythia.readString("Zprime:coup2gen4 = on");
  pythia.readString("32:onMode = off");
  pythia.readString("32:onIfAny = 17");
  pythia.readString("17:m0 = 100.");
  pythia.readString("17:mayDecay = off");
  pythia.readString("PhaseSpace:mHatMin = 200.");
  pythia.init();

  // Book histograms.
  Hist mll("invariant mass of tau' pair", 100, 0., 1000.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    // Fill Z' mass; Z' in position 5.
    mll.fill( pythia.event[5].m() );

  // End of event loop. Statistics. Hostogram. Done.
  }
  pythia.stat();
  cout << mll;
  return 0;
}
