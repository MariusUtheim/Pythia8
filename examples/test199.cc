// test199.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test program of elastic scattering in various frameworks.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate.
  int nEvent = 100000;
  int nAbort = 5;

  // Book histograms.
  Hist tElLog[5];
  tElLog[0].book("elastic log10(|t|) spectrum, mode = 0", 80, -6., 2.);
  tElLog[1].book("elastic log10(|t|) spectrum, mode = 1", 80, -6., 2.);
  tElLog[2].book("elastic log10(|t|) spectrum, mode = 2", 80, -6., 2.);
  tElLog[3].book("elastic log10(|t|) spectrum, mode = 3", 80, -6., 2.);
  tElLog[4].book("elastic log10(|t|) spectrum, mode = 4", 80, -6., 2.);
  Hist   tElLogS("elastic log10(|t|) spectrum, mode sum", 80, -6., 2.);
  double time[5];

  // Loop through five elastic modes.
  for (int iMode = 0; iMode < 5; ++iMode) {
    clock_t start = clock();

    // Generator.
    Pythia pythia;

    // Set up run.
    pythia.readString("Beams:idA = -2212");
    pythia.readString("SoftQCD:elastic = on");
    pythia.settings.mode("SigmaTotal:mode", iMode);
    pythia.readString("SigmaElastic:Coulomb = on");
    pythia.readString("SigmaElastic:tAbsMin = 4e-5");
    pythia.readString("Next:numberCount = 100000");
    pythia.readString("PartonLevel:all = off");

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

      // Study t distribution of elastic events.
      double tAbsL = log10(abs(pythia.info.tHat()));
      tElLog[iMode].fill( tAbsL);

    // End of event loop.
    }

    // Final statistics and normalize histograms.
    pythia.stat();
    tElLog[iMode] *= pythia.info.sigmaGen() * 10. / nEvent;

    // Time usage.
    clock_t stop = clock();
    time[iMode] = double(stop - start) / double(CLOCKS_PER_SEC);
  }
  for (int iMode = 0; iMode < 5; ++iMode) cout << " mode " << iMode
     << " took " << fixed << setprecision(3) << time[iMode] << " s" << endl;

  // Print t spectra before and after normalization to average.
  for (int iMode = 0; iMode < 5; ++iMode) cout << tElLog[iMode];
  for (int iMode = 0; iMode < 5; ++iMode) tElLogS += 0.2 * tElLog[iMode];
  for (int iMode = 0; iMode < 5; ++iMode) tElLog[iMode] /= tElLogS;
  for (int iMode = 0; iMode < 5; ++iMode) cout << tElLog[iMode];

  // Done.
  return 0;
}
