// test195.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check issues in K_S eta spcrum brought up by Deepak Kar.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of event. Process: minbias or ttbar.
  int nEvent = 100000;
  bool doTop = true;

  // Generator. Process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  if (doTop) {
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
  } else {
    pythia.readString("SoftQCD:nonDiffractive");
  }
  pythia.readString("Next:numberShowEvent = 0");

  // Restore earlier behaviour for check.
  pythia.readString("TimeShower:recoilDeadCone = off");

  // Histogram.
  Hist mulC("charged multiplicity", 100, -0.5, 799.5);
  Hist mulK("charged multiplicity", 100, -0.5, 99.5);
  Hist etaC("charged pseudorapidity", 25, 0., 5.);
  Hist etaK("K_S^0 pseudorapidity",   25, 0., 5.);

  // Begin event loop. Generate event. Skip if error.
  pythia.init();
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze event and fill histogram.
    int nCharged = 0;
    int nKaon = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
        ++nCharged;
        etaC.fill( abs(pythia.event[i].eta()) );
      }
      if (pythia.event[i].id() == 310) {
        ++nKaon;
        etaK.fill( abs(pythia.event[i].eta()) );
      }
    }
    mulC.fill( nCharged );
    mulK.fill( nKaon );


  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  mulC /= 4. * nEvent;
  mulK /= nEvent;
  etaC /= 0.2 * nEvent;
  etaK /= 0.2 * nEvent;
  cout << mulC << etaC << mulK << etaK;
  return 0;
}
