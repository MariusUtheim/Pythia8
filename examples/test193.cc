// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readFile("test193.cmnd");

  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  int nEvent = pythia.settings.mode("Main:numberOfEvents");

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    }
    mult.fill( nCharged );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  mult           *= pythia.info.sigmaGen() / pythia.info.nAccepted();

 cout << mult ;

  return 0;
}
