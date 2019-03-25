// main75.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program to study jets in Dark Matter production.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator. Process selection. Initialization. Event shorthand.
  Pythia pythia;
  pythia.readFile("main75.cmnd");
  pythia.init();
  Event& process = pythia.process;


  // Histograms. Error counter.
  Hist pTj("dN/dpTj", 100, 0., 100.);
  Hist mRec("mRec", 100, 0., 1000.);
  int iErr = 0;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) {
      if (++iErr < 100) continue;
      else {
        cout << "Too many errors" << endl;
        break;
      }
    }

    // Invariant mass of DM system.
    Vec4 mRes = process[5].p() + process[6].p();
    mRec.fill(mRes.mCalc());

    // Keep track of missing ET.
    Vec4 missingETvec;

    // Loop over event record to decide what to pass to FastJet.
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final state only.
      if (!pythia.event[i].isFinal()) continue;

      // No neutrinos or DM.
      if ( pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14
        || pythia.event[i].idAbs() == 16 || pythia.event[i].idAbs() == 52)
        continue;

      // Only |eta| < 3.6.
      if (abs(pythia.event[i].eta()) > 3.6) continue;

      // Missing ET.
      missingETvec += pythia.event[i].p();
    }


  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << pTj;

  // Done.
  return 0;
}
