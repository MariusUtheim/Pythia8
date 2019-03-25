// test166.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// gamma + gamma -> H -> gamma + gamma with unresolved p's.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Energy and process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HiggsSM:gmgm2H = on");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 22 22");

  // Force heavy Higgs to be narrow.
  pythia.readString("25:m0 = 750.");
  pythia.readString("25:mWidth = 1.");
  pythia.readString("25:doForceWidth = on");

  // No event printout.
  pythia.readString("Next:numberShowEvent = 5");

  // Set up photon beams with pointlike flux.
  PDF* pdfAPtr = new ProtonPoint ( 2212);
  PDF* pdfBPtr = new ProtonPoint ( 2212);
  pythia.setPDFPtr ( pdfAPtr , 0 );
  pythia.readString("PartonLevel:MPI = off");
  // Do not allow FSR if you want to avoid gamma -> f fbar
  // in Higgs decay.
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:ISR = off");
  // The kT spectrum you get is not the right one, so misleading.
  //pythia.readString("BeamRemnants:primordialKT = off");

  // Switch to force either or both remnants to be original proton.
  pythia.settings.addMode("BeamRemnants:unresolvedHadron",
    0, true, true, 0, 3);
  // Force unresolved: 0 = none, 1 = side 1, 2 = side 2, 3 = both.
  pythia.readString("BeamRemnants:unresolvedHadron = 1");

  // Initialize.
  pythia.init();

  // Histogram.
  Hist massH( "mass of Higgs", 100, 0., 1000.);
  Hist pTH( "pT of Higgs", 100, 0., 100.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze event.
    massH.fill( pythia.info.mHat() );
    int iH = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 25) iH = i;
    pTH.fill( pythia.event[iH].pT() );

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << massH << pTH;

  // Done.
  return 0;
}
