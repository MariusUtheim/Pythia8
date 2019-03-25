// test205.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// pT spectrum of primary gluons in processes.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("SoftQCD:nondiffractive = on");
  pythia.readString("BeamRemnants:primordialKT = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("MultipartonInteractions:alphaSorder = 0");
  pythia.readString("HadronLevel:all = off");
  pythia.init();
  Hist pTa("pT primary quarks&gluons", 100, 0., 10.);
  Hist pTg("pT primary gluons", 100, 0., 10.);
  Hist pT23("pT primary gluons status 23", 100, 0., 10.);
  Hist pT33("pT primary gluons status 33", 100, 0., 10.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find all primary gluons and fill their pT.
    int id, status;
    double pT;
    for (int i = 0; i < pythia.event.size(); ++i) {
      id     = pythia.event[i].id();
      status = pythia.event[i].statusAbs();
      pT     = pythia.event[i].pT();
      if (status == 23 || status == 33) pTa.fill( pT );
      if (id == 21 && (status == 23 || status == 33)) pTg.fill( pT );
      if (id == 21 && status == 23) pT23.fill( pT );
      if (id == 21 && status == 33) pT33.fill( pT );
    }

  // End event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << pTa << pTg << pT23 << pT33;
  return 0;
}
