// test244.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Simple illustration V-A-like decay of fictitious new particle.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Set number of events to generate and to list.
  int nEvent = 10000;
  int nList = 5;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event = pythia.event;

  // Data for new particle 50.
  pythia.readString("50:new = N2 N2 2 0 0 10.0 0.0 0.0 0.0 10.0 0 1 0 1 0");
  pythia.readString("50:addChannel = 1 0.50 94  1 -2 -13");
  pythia.readString("50:addChannel = 1 0.50 94 -1  2  13");
  pythia.readString("50:mayDecay = on");

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Optionally switch off resonance decays, or only showers in them.
  //pythia.readString("ProcessLevel:resonanceDecays = off");
  //pythia.readString("PartonLevel:FSRinResonances = off");

  // Optionally switch off ordinary decays.
  //pythia.readString("HadronLevel:Decay = off");
  pythia.readString("111:mayDecay = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Histogram.
  Hist emu( "muon energy spectrum", 100, 0., 5.001);
  Hist nHad( "number of primary hadrons", 20, -0.5, 19.5); 

  // Initialize.
  pythia.init();

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Reset event record to allow for new event.
    event.reset();

    // Store the particle in the event record.
    event.append( 50, 1, 0, 0, 0., 0., 0., 10., 10.);

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // List first few events. Histogram info.
    if (iEvent < nList) event.list();
    emu.fill( event[2].e() );
    int nh = 0;
    for (int i = 0; i < event.size(); ++i) 
      if (event[i].statusAbs()/10 == 8) ++nh;
    nHad.fill( nh );

  // End of event loop.
  }

  // Print statistics, histograms and done.
  pythia.stat();
  cout << emu << nHad;

  // Done.
  return 0;
}
