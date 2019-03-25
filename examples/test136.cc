// test136.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Excited fermion multibody decays, for Simone Amoroso.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("ExcitedFermion:dg2dStar = on");
  pythia.readString("ExcitedFermion:Lambda = 3800.");
  pythia.readString("4000001:m0=3800");
  pythia.readString("4000001:oneChannel = 1 1 101  22 2 -2 1");
  //pythia.readString("4000001:oneChannel = 1 1 101  22 1 21");
  pythia.readString("Check:nErrList = 1");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("HadronLevel:all = off");
  //pythia.readString("BeamRemnants:primordialKT = off");
  pythia.init();
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  return 0;
}
