// test176.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check new generation of exponential pT spectrum.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Choose parameters of q qbar system.
  int nEvent = 100000;
  int id = 2;
  double ee = 50.0;
  double mm = 0.;
  double pp = sqrt( ee * ee - mm * mm);

  // Generator; shorthand for event.
  Pythia pythia;
  Event& event      = pythia.event;

  // Switch off ProcessLevel and decays.
  pythia.readString("ProcessLevel:all = off");
  pythia.readString("HadronLevel:Decay = off");

  // Exponential decay spectrum.
  pythia.readString("StringPT:thermalModel = on");

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist pTnormal("pT spectrum normal primary hadrons", 100, 0., 2.5);
  Hist pTendpt("pT spectrum endpoint primary hadrons", 100, 0., 2.5);

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Set up  a q qbar system to be hadronized.
    event.reset();
    event.append(  id, 23, 101,   0, 0., 0.,  pp, ee, mm);
    event.append( -id, 23,   0, 101, 0., 0., -pp, ee, mm);

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Fill pT spectra of primary hadrons.
    pTendpt.fill( event[3].pT() );
    for (int i = 4; i < event.size() - 1; ++i) pTnormal.fill ( event[i].pT() );
    pTendpt.fill( event[event.size() - 1].pT() );

  // End of event loop.
  }

  // Print statistics, histograms and done.
  pythia.stat();
  cout << pTnormal << pTendpt;

  // Done.
  return 0;
}
