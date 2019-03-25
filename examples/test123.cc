// test123.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Inelastisity distribution in neutrino interactions.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 10000;

  // Loop over nu_mu and nu_mubar beams.
  for (int iBeam = 0; iBeam < 2; ++iBeam) {

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Beam parameters.
  if (iBeam == 0) pythia.readString("Beams:idA =  14");
  else            pythia.readString("Beams:idA = -14");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:frameType = 3");
  pythia.readString("Beams:pzA = 100000.");
  pythia.readString("Beams:pzB = 0.");

  // Process selection.
  pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");

  // Switch off ISR, FSR, primordial kT: not yet implemented.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("Check:nErrList = 5");

  // Change process PDF or factorization scale.
  //pythia.readString("PDF:pSet = 15");
  //pythia.readString("SigmaProcess:factorScale2 = 5");

  // Initialization.
  pythia.init();

  // Histogram.
  Hist inel("inelasticity", 20, 0., 1.);
  Hist pTmu("pT of muon", 100, 0., 200.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    //cout << " now begin event no " << iEvent << endl;
    if (!pythia.next()) continue;

    // Calculate and histogram inelasticity.
    double x = 1. - event[5].e() / event[1].e();
    inel.fill( x );
    pTmu.fill( event[5].pT() );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  double sigma = pythia.info.sigmaGen();
  inel *= sigma * 1e7 * 20. / nEvent;
  pTmu *= sigma * 1e7 * 2. / nEvent;
  cout << inel << pTmu;

  // End loop over nu_mu and nu_mubar beams.
  }

  // Done.
  return 0;
}
