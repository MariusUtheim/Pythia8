// test165.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Read in commands.
  pythia.readString("Beams:idA = 22");
  pythia.readString("Beams:idB = 22");
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("HardQCD:all = on");
  pythia.readString(" PhaseSpace:pTHatMin = 10.");
  //pythia.readString("");

  // Settings to be used in the main program.
  int nEvent = 1000;
  int nAbort = 10;

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist pTHat("pThat of hard interaction", 100, 0., 50.);
  Hist nChg(  "charged multiplicity", 100, -1., 199.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Study event.
    double pTH = pythia.info.pTHat();
    pTHat.fill( pTH);
    int nCh = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nCh;
    nChg.fill( nCh);

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.stat();
  cout << pTHat << nChg;

  // Done.
  return 0;
}
