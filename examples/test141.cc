// test141.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Initialize several Pythia instances from each other.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  /*
  // Loop to create Pythia instances, for timing purposes.
  for (int iLoop = 0; iLoop < 100; ++iLoop) {
    // Pythia pythiaLoop;
    Pythia pythiaLoop( pythia.settings, pythia.particleData);
    // pythiaLoop.init();
  }
  */

  // LHC initialization. Histogram.
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Copy of first Pythia.
  Pythia pythia2( pythia.settings, pythia.particleData);
  pythia2.init();

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia2.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia2.event.size(); ++i)
      if (pythia2.event[i].isFinal() && pythia2.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

  // End of event loop. Statistics.
  }
  pythia2.stat();

  // Histogram.
  cout << mult;

  // Compare particle data files (use diff manually).
  //pythia.particleData.listFF("out141a");
  //pythia2.particleData.listFF("out141b");

  return 0;
}
