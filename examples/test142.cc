// test142.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check double onium production.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 10000;

  // Generator.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("Next:numberShowProcess = 5");

  // Process selection.
  pythia.readString("Charmonium:gg2doubleccbar(3S1)[3S1(1)] = on,on,on");
  pythia.readString("Charmonium:qqbar2doubleccbar(3S1)[3S1(1)] = on,on,on");
  pythia.readString(
    "Bottomonium:gg2doublebbbar(3S1)[3S1(1)] = on,on,on,on,on,on");
  pythia.readString(
    "Bottomonium:qqbar2doublebbbar(3S1)[3S1(1)] = on,on,on,on,on,on");

  // For cross sections possible to simplify generation.
  pythia.readString("PartonLevel:all = off");
  pythia.readString("HadronLevel:all = off");

  // Initialization.
  pythia.init();

  // Histograms.
  Hist mProc( "invariant mass of hard process", 100, 0., 100.);
  Hist pTProc( "transverse momentum of hard process", 100, 0., 20.);

  // Begin event loop. Generate events. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Study process invariant mass.
    mProc.fill( pythia.info.mHat() );
    pTProc.fill( pythia.info.pTHat() );

  // End of event loop. Statistics. Done.
  }
  pythia.stat();
  cout << mProc << pTProc;
  return 0;
}
