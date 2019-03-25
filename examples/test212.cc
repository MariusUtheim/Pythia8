// test212.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.


#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Number of events.
  int nEvent = 10000;

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("MultipartonInteractions:processLevel = 2");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 399.5);
  Hist dndyP("dn/dy charged positive", 100, 0., 10.);
  Hist dndyN("dn/dy charged negative", 100, 0., 10.);
  Hist dndyD("dn/dy charged difference", 100, 0., 10.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
      ++nCharged;
      double yNow = pythia.event[i].y();
      if (yNow > 0.) dndyP.fill(yNow);
      if (yNow < 0.) dndyN.fill(-yNow);
    }
    mult.fill( nCharged );

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  dndyP *= 10. / nEvent;
  dndyN *= 10. / nEvent;
  dndyD = dndyP - dndyN;
  cout << mult << dndyP << dndyN << dndyD;
  return 0;
}
