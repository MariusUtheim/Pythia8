// test114.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Deeply Inelastic Scattering - preliminary.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Set up for HERA incoming beams.
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:frameType = 2");
  pythia.readString("Beams:eA = 820.");
  pythia.readString("Beams:eB = 27.5");

  // No QED radiation off e (data usually corrected to this).
  pythia.readString("PDF:lepton = off");

  // No ISR, FSR. MPI and primordial kT not present so need no switch.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("Check:nErrList = 3");

  // DIS processes, both neutral and charged current.
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");

  // Kinematics cut; not quite suited for DIS thinking.
  pythia.readString("PhaseSpace:pTHatMin = 10.");

  // Initialize.
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    //cout << "\n Begin event number " << iEvent << endl;
    if (!pythia.next()) continue;

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
