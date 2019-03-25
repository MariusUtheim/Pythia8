// test157.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Read in commands from external file.
  pythia.readString("SoftQCD:elastic = on");

  // Extract settings to be used in the main program.
  int    nEvent = 100000;
  int    nAbort = 10;

  // Initialize.
  pythia.init();

  // Book histograms: elastic/diffractive.
  Hist tSpecEl("elastic |t| spectrum",              100, 0., 1.);
  Hist tWtEl(  "elastic |t| spectrum, with weight", 100, 0., 1.);
  Hist tSpecElLog("elastic log10(|t|) spectrum",    100, -5., 0.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Study t distribution of elastic events.
    double tAbs = abs(pythia.info.tHat());
    double wt   = pythia.info.weight();
    tSpecEl.fill(tAbs);
    tWtEl.fill(tAbs, wt);
    tSpecElLog.fill(log10(tAbs));

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.stat();
  cout << tSpecEl << tWtEl << tSpecElLog;

  // Done.
  return 0;
}
