// test101.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  //pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("SoftQCD:nonDiffractive = on");
  //pythia.readString("SoftQCD:singleDiffractive = on");
  pythia.readString("SigmaTotal:zeroAXB = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("HadronLevel:all = off");
  //pythia.readString("PDF:pSet = 2");
  pythia.readString("Check:nErrList = 2");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  Hist multMPI("number of MPIs", 50, -0.5, 49.5);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );
    // Check if MPI.
    int nMPI = 1;
    for (int i = 1; i < pythia.event.size(); ++i)
      if ( pythia.event[i].status() == -31
        && pythia.event[i - 1].status() != -31) ++nMPI;
    multMPI.fill( nMPI );

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << mult << multMPI;

  // Check output method.
  cout << pythia.settings.output("PartonLevel:ISR")
       << pythia.settings.output("Check:nErrList")
       << pythia.settings.output("Beams:eCM")
       << pythia.settings.output("PDF:pSet")
       << pythia.settings.output("Charmonium:gg2ccbar(3PJ)[3S1(8)]g")
       << pythia.settings.output("Charmonium:states(3PJ)")
       << pythia.settings.output("Charmonium:O(3PJ)[3S1(8)]")
       << pythia.settings.output("junkTag");
  pythia.readString("Beams:eCM = ?");
  pythia.readString("SoftQCD:nonDiffractive = ?");
  pythia.readString("Check:nErrList = ?");

  return 0;
}
