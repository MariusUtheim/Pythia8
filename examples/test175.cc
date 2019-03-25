// test175.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program for (anti)neutron PDF access.
// To access LHAPDF6 via the LHAPDF5 Fortran interface use
// ./configure --with-lhapdf6=/Users/torbjorn/code/lhapdf
// or equivalent on your installation.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Beam combinations: 1 = pp, 2 = pn, 3 = nn,
  // 4 = ppbar, 5 = pnbar, 6 = nnbar, 7 = nbarnbar.
  double sigma[8];
  string label[8] = { " ", "pp", "pn", "nn", "ppbar", "pnbar", "nnbar",
    "nbarnbar" };
  for (int combi = 1; combi < 8; ++combi) {

    // Generator and beam combination.
    Pythia pythia;
    if (combi < 6  && combi != 3) pythia.readString("Beams:idA =  2212");
    if (combi == 3 || combi == 6) pythia.readString("Beams:idA =  2112");
    if (combi == 7)               pythia.readString("Beams:idA = -2112");
    if (combi == 1)               pythia.readString("Beams:idB =  2212");
    if (combi == 2 || combi == 3) pythia.readString("Beams:idB =  2112");
    if (combi == 4)               pythia.readString("Beams:idB = -2212");
    if (combi > 4)                pythia.readString("Beams:idB = -2112");

    // Set up pure gamma* production => sensitive to e_q^2,
    // and ensure large x values so valence quarks dominate.
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("WeakZ0:gmZmode = 1");
    pythia.readString("Beams:eCM = 200.");
    pythia.readString("PhaseSpace:mHatMin = 40.");
    pythia.readString("PhaseSpace:mHatMax = 60.");

    // Only generate hard process. Select LHAPDF6 set.
    pythia.readString("PartonLevel:all = off");
    pythia.readString("PDF:pSet = LHAPDF6:MRST2007lomod.LHgrid");

    // Initialize, generate events and print statistics.
    pythia.init();
    for (int iEvent = 0; iEvent < 10000; ++iEvent) pythia.next();
    pythia.stat();
    sigma[combi] = pythia.info.sigmaGen();
  }

  // Print summary table.
  cout << "\n\n combi   label     sigmaGen" << endl;
  for (int combi = 1; combi < 8; ++ combi) cout << setw(4) << combi
    << setw(10) << label[combi] << setw(14) << scientific
    << setprecision(4) << sigma[combi] << endl;

  return 0;
}
