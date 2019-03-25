// test177.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study consequences of new Q2min cut, implying new z = cos(thetaHat)
// implementation.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  //pythia.readString("HardQCD:qg2qg = on");
  //pythia.readString("HardQCD:qq2qq = on");
  pythia.readString("HardQCD:gg2gg = on");
  //pythia.readString("PhaseSpace:mHatMin = 100.");
  //pythia.readString("PhaseSpace:mHatMax = 120.");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  //pythia.readString("PhaseSpace:pTHatMax = 25.");
  //pythia.readString("PhaseSpace:Q2Min = 400.");
  pythia.readString("PartonLevel:all = off");
  pythia.init();

  // Histograms.
  Hist pTH("pT", 100, 0., 100.);
  Hist QH("Q = sqrt(-tHat)", 100, 0., 100.);
  Hist ctH("cos(thetaHat)", 100, -1., 1.);
  Hist mH("mHat", 100, 0., 1000.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100000; ++iEvent) {
    if (!pythia.next()) continue;

    // Fill histograms.
    double pT = pythia.info.pTHat();
    double sH = pythia.info.sHat();
    double tH = pythia.info.tHat();
    double uH = pythia.info.uHat();
    pTH.fill( pT );
    QH.fill( sqrt(abs(tH)) );
    ctH.fill( (tH - uH) / sH );
    mH.fill( sqrt(sH) );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << pTH << QH << ctH << mH;

  // Done.
  return 0;
}
