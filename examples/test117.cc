// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch, based on main01.cc.
// This program illustrates how HepMC can be interfaced to Pythia8.
// It studies the charged multiplicity distribution at the LHC.
// HepMC events are output to the hepmcout41.dat file.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:frameType = 2");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eA = 100.");
  pythia.readString("Beams:eB = 100.");

  pythia.readString("Init:showChangedSettings = on");

  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
  pythia.readString("Zprime:gmZmode = 3");

  pythia.readString("32:m0 = 4.");
  pythia.readString("32:mMin = 0.1");
  pythia.readString("32:mMax = 1000.");
  pythia.readString("PhaseSpace:mHatMin = 0.1");
  pythia.init();

  Hist mZp( "mass of Zprime", 100, 0., 10.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
    int iZp = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 32) iZp = i;
    }
    mZp.fill( pythia.event[iZp].m() );

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << mZp;

  // Done.
  return 0;
}
