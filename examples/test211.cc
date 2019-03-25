// test211.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Impact-parameter setting for diffraction.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate and abort criterium.
  int nEvent = 100;
  int nAbort = 3;

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up generation.
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 8000.");
  //pythia.readString("SoftQCD:all = on");
  pythia.readString("SoftQCD:singleDiffractive = on");
  pythia.readString("SoftQCD:doubleDiffractive = on");
  //pythia.readString("SoftQCD:centralDiffractive = on");
  //pythia.readString("SoftQCD:nonDiffractive = on");
  //pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize.
  pythia.init();

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Extract event classification.
    int code = pythia.info.code();
    double bMPI = pythia.info.bMPI();
    double enhance = pythia.info.enhanceMPI();
    int nMPI = pythia.info.nMPI();
    double m1 = event[3].m();
    double m2 = event[4].m();
    cout << scientific << setprecision(3) << setw(10) << bMPI
         << "   " << setw(3) << code << "   " << setw(10) << m1
         << "   " << setw(10) << m2 << "   " << setw(3) << nMPI
         << "   " << setw(10) << enhance << endl;

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.stat();

  // Done.
  return 0;
}
