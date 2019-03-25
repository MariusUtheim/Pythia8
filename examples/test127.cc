// test127.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug production of high-pT particles from junction bugs.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

// Histogram made global.
Hist rootval("log10(square root argument all pairs)", 100, -10., 10.);
Hist rootval2("log10(square root argument selected)", 100, -10., 10.);

int main() {

  // Number of events.
  int nEvent = 100000;

  // Generator. Process selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("SoftQCD:inelastic = on");
  //pythia.readString("ColourReconnection:mode = 1");  // More junctions!
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << rootval << rootval2;
  return 0;
}

// in StringFragmentation.cc, after #include, add

//extern Pythia8::Hist rootval;
//extern Pythia8::Hist rootval2;

// at line 971, just before M2MINJRF check, add

//    rootval.fill( log10((pWTinJRF[0] + pWTinJRF[1]).m2Calc()) );
//    rootval.fill( log10((pWTinJRF[0] + pWTinJRF[2]).m2Calc()) );
//    rootval.fill( log10((pWTinJRF[1] + pWTinJRF[2]).m2Calc()) );

// at line 1297, just before M2MINJRF check, add

//      rootval2.fill( log10(m2j + m2k + 2. * pjpk) );
