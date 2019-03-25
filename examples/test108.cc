// tesy108.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Generate a predetermined second hard interaction.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Loop over pT0 to study dependence.
  for (int iPT = 0; iPT < 4; ++iPT) {

  // Generator.
  Pythia pythia;

  // Select first hard process (just a small sample of possibilities).
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("Top:all = on");
  //pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("WeakSingleBoson:ffbar2W = on");

  // Select second hard process (complete list of options).
  pythia.readString("SecondHard:generate = on");
  //pythia.readString("SecondHard:TwoJets = on");
  //pythia.readString("SecondHard:PhotonAndJet = on");
  //pythia.readString("SecondHard:TwoPhotons = on");
  //pythia.readString("SecondHard:SingleGmZ = on");
  pythia.readString("SecondHard:SingleW = on");
  //pythia.readString("SecondHard:TwoBJets = on");

  // Kinematics cuts, common for the two.
  pythia.readString("PhaseSpace:mHatMin = 40.");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  // Remove uninteresting aspects for sigma_eff.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("HadronLevel:all = off");

  // Vary pTmin to figure out dependence.
  if (iPT == 0) pythia.readString("MultipartonInteractions:pT0Ref = 1.2");
  if (iPT == 1) pythia.readString("MultipartonInteractions:pT0Ref = 1.6");
  if (iPT == 2) pythia.readString("MultipartonInteractions:pT0Ref = 2.0");
  if (iPT == 3) pythia.readString("MultipartonInteractions:pT0Ref = 2.4");

  // Initialize for LHC at 8 TeV.
  pythia.readString("Beams:eCM = 8000.");
  pythia.init();

  // Histogram.
  Hist nMult("number of multiparton interactions", 100, -0.5, 99.5);
  Hist bMore("b enhancement factor",    100, 0., 10.);

  // Generate events.
  for (int iev = 0; iev < 1000; ++iev) {
    pythia.next();

    // Histogram multiparton interactions
    double nMPI = pythia.info.nMPI();
    nMult.fill( nMPI );
    bMore.fill( pythia.info.enhanceMPI() );

  }

  // Compare full statistics listing with what is set in info.
  pythia.stat();
  cout << scientific << setprecision(3) << "\n From pythia.info: sigma = "
       << pythia.info.sigmaGen() << " +- " << pythia.info.sigmaErr()
       << endl;

  // Print histograms.
  cout << nMult << bMore;
  }

  // Done.
  return 0;
}
