// test118.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test jet spectra in diffraction.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//============================================================================

int main() {

  // Number of events.
  int nEvent = 10000000;

  // Generator. Tevatron QCD jets.
  Pythia pythia;
  pythia.readString("Beams:idB = -2212");
  pythia.readString("Beams:eCM = 1800.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 8.");

  // Allow hard diffraction.
  pythia.readString("Diffraction:doHard = on");
  pythia.readString("Diffraction:sampleType = 1");

  // Switch off irrelevant parts of the generation.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Next:numberCount = 0");

  // Initialize. Book histograms.
  pythia.init();
  Hist pTallH("pT spectrum all events",         100, 0., 100.);
  Hist pTdifH("pT spectrum diffractive events", 100, 0., 100.);
  Hist pTratH("pT fraction diffractive events", 100, 0., 100.);
  Hist xAllH( "all x values",                   100, 0., 0.2);
  Hist xDifH( "diffractive x values",           100, 0., 0.2);
  Hist xRatH( "fraction diffractive x values",  100, 0., 0.2);
  Hist xRecH( "recoiler x values",              100, 0., 0.2);
  Hist xReRH( "fraction recoiler x values",     100, 0., 0.2);
  Hist xPomH( "Pomeron x values",               100, 0., 1.0);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find diffractive character.
    bool isHDA = pythia.info.isHardDiffractiveA();
    bool isHDB = pythia.info.isHardDiffractiveB();
    bool isHD  = isHDA || isHDB;

    // Histogram pT and x of all hard interactions.
    double pTHat = pythia.info.pTHat();
    pTallH.fill( pTHat );
    double x1 = pythia.info.x1();
    double x2 = pythia.info.x2();
    xAllH.fill( x1 );
    xAllH.fill( x2 );

    // Histogram pT, x and xPom values of diffractive events.
    if (isHD) {
      double xD   = isHDB ? x1 : x2;
      double xR   = isHDB ? x2 : x1;
      double xPom = isHDB ? pythia.info.xPomeronA() : pythia.info.xPomeronB();
      xPomH.fill( xPom );
      if (xPom > 0.035 && xPom < 0.095) {
        pTdifH.fill( pTHat );
        xDifH.fill( xD );
        xRecH.fill( xR );
      }
    }

  // End of event loop. Statistics. Histograms. Done.
  }
  pythia.stat();
  pTratH = pTdifH / pTallH;
  xRatH  = xDifH / xAllH;
  xReRH  = xRecH / xAllH;
  cout << pTallH << pTdifH << pTratH << xAllH << xDifH << xRatH
       << xRecH << xReRH << xPomH;
  return 0;
}
