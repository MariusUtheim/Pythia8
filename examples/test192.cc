// test192.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Compare charged and MPI mutiplicity for nondiffractive.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("SoftQCD:nondiffractive = on");
  pythia.init();
  Hist multChg("charged multiplicity", 100, -0.5, 799.5);
  Hist multMPI("MPI multiplicity", 100, -0.5, 99.5);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nChg = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nChg;
    multChg.fill( nChg );
    multMPI.fill( pythia.info.nMPI() );

  // End of event loop. Statistics. Histograms. Done.
  }
  pythia.stat();
  cout << multChg << multMPI;
  return 0;
}
