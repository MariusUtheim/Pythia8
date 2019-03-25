// test129.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug four-junction event (Gabriel Magill).

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;

  // Set up for Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  //pythia.readString("Beams:LHEF = test1294jun.lhe");
  //pythia.readString("Beams:LHEF = test129long.lhe");
  pythia.readString("Beams:LHEF = monotop_part.lhe");

  // Optional commands.
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("HadronLevel:all = off");

  // Initialize.
  pythia.init();

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; ; ++iEvent) {

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    // Optional checks.
    if (iEvent < 30) pythia.event.listJunctions();

  // End of event loop.
  }

  // Give statistics.
  pythia.stat();

  // Done.
  return 0;
}
