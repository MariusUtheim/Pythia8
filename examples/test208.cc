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
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Fragmentation:setVertices = on");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
    mult.fill( nCharged );

    /*
    // Search for octet states and check whether hadrons.
    for (int i = 0; i < event.size(); ++i) {
      bool print = false;
      if (event[i].idAbs() > 9900000) print = true;
      if (event[i].isHadron() && (event[i].col() != 0 || event[i].acol() != 0))
        print = true;
      if (print) cout << " code " << event[i].id() << " is called "
         << event[i].name() << " and isHadron = " << event[i].isHadron()
	 << endl;
    }
    */

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << mult;
  return 0;
}
