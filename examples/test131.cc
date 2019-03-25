// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("RHadrons:allow = on");
  pythia.settings.forceParm("RHadrons:maxWidth", 3.);
  pythia.readString("RHadrons:idStop = 6");
  pythia.init();
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {
    if (!pythia.next()) continue;
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  return 0;
}
