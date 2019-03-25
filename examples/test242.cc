// test241.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Compare partonic and hadronic multiplicities

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 100000;
  int nAbort  = 5;

  // Book histograms.
  Hist nChH( "n_charged", 100, -0.5, 99.5);
  Hist nMPIH( "n_MPI", 100, -0.5, 99.5);
  Hist nChMPIH( "n_charged(n_MPI)", 100, -0.5, 99.5);
  Hist nMPIChH( "n_MPI(n_charged)", 100, -0.5, 99.5);

  // Generator setup.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("SoftQCD:nondiffractive = on");
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

    // Count charged hadrons in |eta| < 1. Find number of MPIs.
    int nCh = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()
      && event[i].isCharged() && abs(event[i].eta()) < 1.) ++nCh;
    int nMPI = pythia.info.nMPI();

    // Fill histograms.
    nChH.fill( nCh);
    nMPIH.fill( nMPI);
    nChMPIH.fill( nMPI, nCh);
    nMPIChH.fill( nCh, nMPI);


  // End of event loop. Print histograms.
  }
  pythia.stat();
  nChMPIH /= nMPIH;
  nMPIChH /= nChH;
  cout << nChH << nMPIH << nChMPIH << nMPIChH;

  // Done.
  return 0;
}
