// main37.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File version 3.0 input can be used
// to mix events according to several different event weights.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator
  Pythia pythia;

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = ttbar.lhe");
  //pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
  //pythia.readString("Beams:LHEF = w+_production_lhc_2.lhe");
  pythia.readString("5:m0 = 10.");
  pythia.readString("11:m0 = 10.");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Next:numberShowLHA = 10");
  pythia.readString("Next:numberShowProcess = 10");
  pythia.init();

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      else continue;
    }

  // End of event loop.
  }

  // Give statistics.
  pythia.stat();

  // Done.
  return 0;
}
