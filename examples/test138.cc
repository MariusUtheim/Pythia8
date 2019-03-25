// main25.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// Here the very few events are generated with MadGraph, and illustrate
// more complicated colour topologies.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator
  Pythia pythia;

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:Hadronize = on");

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = test138.lhe");
  //pythia.particleData.listAll();
  //pythia.particleData.listChanged();
  pythia.init();


  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; ; ++iEvent) {
    cout << endl << "Begin event # " << iEvent << endl;

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    //pythia.process.list();
    pythia.event.list(true,true);

  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();

  // Done.
  return 0;
}
