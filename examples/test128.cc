// test128.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study B_s - B_sbar oscillations.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events to generate.
  int nEvent = 50000;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // List first three events.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 3");

  // Extract B0S properties.
  int    idB0S   = 531;
  double mB0S    = pdt.m0(idB0S);
  double tau0B0S = pdt.tau0(idB0S);

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist tauU("decay lifetime unmixed", 100, 0., 2.5);
  Hist tauM("decay lifetime mixed  ", 100, 0., 2.5);

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Set up single B0S.
    event.reset();
    event.append( 531, 1, 0, 0, 0., 0., 0., mB0S, mB0S);
    double tau = tau0B0S * pythia.rndm.exp();
    event[1].tau( tau);

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Find and histogram lifetime.
    bool noOsc = (event[2].statusAbs() == 91);
    if (noOsc) tauU.fill( tau / tau0B0S);
     else      tauM.fill( tau / tau0B0S);

  // End of event loop.
  }

  // Print statistics, histograms and done.
  pythia.stat();
  cout << tauU << tauM;

  // Done.
  return 0;
}
