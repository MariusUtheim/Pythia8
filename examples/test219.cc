// test219.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
  pythia.readString("32:m0 = 450.");
  pythia.readString("PhaseSpace:mHatMin = 400.");
  pythia.readString("PhaseSpace:mHatMax = 500.");
  pythia.readString("32:tau0 = 10.");
  pythia.init();
  Hist pTZp("dN/dpTZp", 100, 0., 1000.);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent == 0) pythia.event.list(true);
    // Loop over particles in event. Find last Z0p copy. Fill its pT.
    int iZp = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 32) iZp = i;
    pTZp.fill( pythia.event[iZp].pT() );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << pTZp;
  return 0;
}
