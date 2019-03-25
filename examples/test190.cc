// test190.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies secondary c and b production in light-quark LEP1 events.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator and event record.
  Pythia pythia;
  Event& event = pythia.event;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3");

  // Do not need to hadronize for these studies.
  pythia.readString("HadronLevel:all = off");

  // Test option with ME corections switched off after first.
  pythia.readString("TimeShower:MEafterFirst = off");

  // LEP1 initialization at Z0 mass.
  pythia.readString("Beams:idA =  11");
  pythia.readString("Beams:idB = -11");
  double mZ = pythia.particleData.m0(23);
  pythia.settings.parm("Beams:eCM", mZ);
  pythia.init();

  // Statistics on c, b and gproduction.
  int nC = 0;
  int nB = 0;
  int nG = 0;

  // Begin event loop. Generate event. Skip if error. List first few.
  for (int iEvent = 0; iEvent < 1000000; ++iEvent) {
    if (!pythia.next()) continue;

    // Search for final c and b quarks.
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].id() ==  4 && event[i].isFinal()) ++nC;
      if (event[i].id() ==  5 && event[i].isFinal()) ++nB;
      if (event[i].id() == 21 && event[i].isFinal()) ++nG;
    }

  // End of event loop. Statistics.
  }
  pythia.stat();
  cout << " number of c quarks = " << nC << endl;
  cout << " number of b quarks = " << nB << endl;
  cout << " number of gluons = " << nG << endl;

  // Done.
  return 0;
}
