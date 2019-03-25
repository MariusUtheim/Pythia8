// test103.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check/debug of new compositeness processes,
// based on code from Olga (Olya) Igonkina, December 2014.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // 1 = new processes; 2 = production via Z' machinery.
  int mode = 1;
  int nEvent = 10000;

  // Generator.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("Next:numberShowProcess = 5");

  // Process selection.
  if (mode == 1) {
    //pythia.readString("ExcitedFermion:qqbar2muStarmu = on");
    pythia.readString("ExcitedFermion:qqbar2muStarmuStar = on");
    pythia.readString("4000013:m0 = 2000.");
    pythia.readString("ExcitedFermion:Lambda = 4000.");
    //pythia.readString("4000013:doForceWidth = on");
    //pythia.readString("4000013:mWidth = 1.");
  } else {
    pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
    pythia.readString("32:onMode = off");
    pythia.readString("32:addChannel = 1 1. 100 4000013 -4000013");
    pythia.readString("32:mMin = 800.");
    pythia.readString("Zprime:gmZmode = 4");
  }

  // For cross sections possible to simplify generation.
  pythia.readString("PartonLevel:all = off");
  pythia.readString("HadronLevel:all = off");

  // Initialization.
  pythia.init();
  //pythia.particleData.list(4000001);
  //pythia.particleData.list(4000013);

  // Histograms.
  Hist mProc( "invariant mass of hard process", 100, 0., 8000.);
  Hist mStar( "mass of excited fermion", 100, 0., 8000.);

  // Begin event loop. Generate events. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Study process invariant mass.
    mProc.fill( pythia.info.mHat() );
    if (mode == 1) mStar.fill( pythia.process[5].m() );
    if (mode == 1) mStar.fill( pythia.process[6].m() );
    if (mode == 2) mStar.fill( pythia.process[6].m() );
    if (mode == 2) mStar.fill( pythia.process[7].m() );

  // End of event loop. Statistics. Done.
  }
  pythia.stat();
  cout << mProc << mStar;
  return 0;
}
