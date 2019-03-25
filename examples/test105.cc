// test105.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Consequences of suppressing primordial kT at large rapidities.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate.
  int nEvent = 100000;

  // Histograms.
  Hist dndy0("dn/d|y| charged par = 0.",         100, 0., 10.);
  Hist dndy1("dn/d|y| charged par = 0.5",        100, 0., 10.);
  Hist dndyr("dn/d|y| charged par = 0.5 / 0.",   100, 0., 10.);
  Hist pTy0( "<pT>(|y|) charged par = 0.",       100, 0., 10.);
  Hist pTy1( "<pT>(|y|) charged par = 0.5",      100, 0., 10.);
  Hist pTyr( "<pT>(|y|) charged par = 0.5 / 0.", 100, 0., 10.);

  // Loop over two suppressions of primordial kT.
  for (int mode = 0; mode < 2; ++mode) {

    // Generator; shorthand for event.
    Pythia pythia;
    Event& event = pythia.event;

    // Select process.
    pythia.readString("SoftQCD:nonDiffractive = on");
    //pythia.readString("SoftQCD:inelastic = on");
    if (mode == 0) pythia.readString("BeamRemnants:reducedKTatHighY = 0.0");
    if (mode == 1) pythia.readString("BeamRemnants:reducedKTatHighY = 0.5");

    // Restrict output.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    // Initialize.
    pythia.init();

    // Begin of event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Skip if failure.
      if (!pythia.next()) continue;

      // Loop to find all final charged particles.
      double yA, pT;
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) {
        yA = abs(event[i].y());
        pT = event[i].pT();
        if (mode == 0) {
          dndy0.fill( yA);
          pTy0.fill( yA, pT);
        } else {
          dndy1.fill( yA);
          pTy1.fill( yA, pT);
        }
      }

    // End of event and mode loops. Statistics.
    }
    pythia.stat();
  }

  // Print histograms.
  pTy0 /= dndy0;
  pTy1 /= dndy1;
  dndy0 *= 10. / nEvent;
  dndy1 *= 10. / nEvent;
  dndyr = dndy1 / dndy0;
  pTyr = pTy1 / pTy0;
  cout << dndy0 << dndy1 << dndyr << pTy0 << pTy1 << pTyr;

  // Done.
  return 0;
}
