// main55.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.


#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 1000;

  // Generator and collision energy.
  Pythia pythia;
  pythia.readFile("test161.cmnd");

  // Initialize. Shorthand for event.
  pythia.init();

  // Histograms.
  Hist mHmot( "original H mass", 100, 250., 300.);
  Hist mHmin( "heavy h mass", 100, 0., 200.);
  Hist mHmax( "light h mass", 100, 0., 200.);


  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

  // Fill histograms.
  double mMax = max( pythia.process[6].m(), pythia.process[7].m() );
  double mMin = min( pythia.process[6].m(), pythia.process[7].m() );
  mHmot.fill( pythia.process[5].m() );
  mHmin.fill( mMax );
  mHmax.fill( mMin );

  // End of event loop. Statistics. Histograms
  }
  pythia.stat();
  cout << mHmot << mHmax << mHmin;

  // Done.
  return 0;
}
