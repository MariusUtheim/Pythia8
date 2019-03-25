// test233.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug SK CR for Marina Beguin.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("test233.cmnd");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();
  int nRecon = 0;
  Hist nRecPart( "number of reconnected partons", 80, -0.5, 79.5);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    //cout << " begin event # = " << iEvent << endl;

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Check if colour reconnection has occured.
    int nRecNow = 0;
    for (int i = 6; i < event.size(); ++i)
    if (event[i].statusAbs() == 79) ++nRecNow;
    nRecPart.fill( nRecNow);
    if (nRecNow > 0) ++nRecon;

  // End of event loop. Statistics.
  }
  pythia.stat();
  cout << nRecPart;
  cout << " Number of reconnected events = " << nRecon << endl;

  // Done.
  return 0;
}
