// test110.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// FSR QED.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 1000;

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 1000.");
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PDF:lepton = off");
  //pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString("TimeShower:allowBeamRecoil = off");
  pythia.readString("23:onMode = off");
  pythia.readString("23:OnIfMatch = 11 -11");
  pythia.init();

  // Histograms.
  Hist nGamHist("number of photons", 100, -0.5, 99.5);
  Hist pTGamHist("relative pT of photons", 100, 0., 2.);

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    double pTHard = pythia.event[5].pT();

    // Find number of final photons and their pTs. Histogram.
    int nGam = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].id() == 22) {
       ++nGam;
       pTGamHist.fill( pythia.event[i].pT() / pTHard );
    }
    nGamHist.fill( nGam );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << nGamHist << pTGamHist;

  // Done.
  return 0;
}
