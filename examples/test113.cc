// test113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Higgs production with subsequent decay to gamma phi.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Energy and process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HiggsSM:all = on");

  // Add gamma + phi decay channel to Higgs and switch off others.
  pythia.readString("25:addchannel = 1 0.00001 101 22 333");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 22 333");

  // Initialize.
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
