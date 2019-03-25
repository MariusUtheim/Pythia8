// test104.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug negative-energy particle.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Set number of events to generate and to list.
  int nEvent = 10000;
  int nList  = 3;
  int nAbort = 10;

  // Properties of event to generate.
  int id    = 1;
  double pp = 0.8;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event      = pythia.event;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Optionally switch off ordinary decays.
  //pythia.readString("HadronLevel:Decay = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Print flawed events.
  //pythia.readString("Check:event = off");
  pythia.readString("Check:nErrList= 5");

  // Initialize.
  pythia.init();

  // Begin of event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Fill event to be hadronized.
    event.reset();
    if (id == 21) {
      event.append( 21, 23, 101, 102, 0., 0.,  pp, pp, 0.);
      event.append( 21, 23, 102, 101, 0., 0., -pp, pp, 0.);
    } else {
      event.append(  id, 23, 101,   0, 0., 0., 0.5 * pp, 0.5 * pp, 0.);
      event.append(  21, 23, 102, 101, 0., 0.,      -pp,       pp, 0.);
      event.append( -id, 23,   0, 102, 0., 0., 0.5 * pp, 0.5 * pp, 0.);
    }

    // Generate events. Quit if failure.
    if (!pythia.next() && ++iAbort > nAbort) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // List first few events.
    if (iEvent < nList) event.list();

    // Loop over all particles to find negative-energy ones.
    for (int i = 0; i < event.size(); ++i) if (event[i].e() < 0.) {
      cout << "\n Error! Particle " << i << " has negative energy!";
      event.list();
    }

  // End of event loop.
  }

  // Print statistics and done.
  pythia.stat();
  return 0;
}
