// test239.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Generate a predetermined second hard interaction.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Key run parameters.
  int nEvent = 10000;
  bool outputLHEF = true;
  bool inputLHEF  = false;
  bool findJets   = true;
  string fileName = "dpsj.lhe";

  // Generator.
  Pythia pythia;
  Event& event   = pythia.event;

  // Create an LHAup object that can access relevant information in pythia.
  // Open a file on which LHEF events should be stored, and write header.
  if (inputLHEF == true) outputLHEF = false;
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  if (outputLHEF) myLHA.openLHEF(fileName);

  // Option with input from LHEF.
  if (inputLHEF) {
    pythia.readString("Beams:frameType = 4");
    pythia.settings.word("Beams:LHEF", fileName);

  // Normal choice with internal generation. CM energy choice.
  } else {
    pythia.readString("Beams:eCM = 8000.");

    // Select first hard process (just a small sample of possibilities).
    pythia.readString("HardQCD:all = on");
    //pythia.readString("Top:all = on");
    //pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    //pythia.readString("WeakSingleBoson:ffbar2W = on");

    // Select second hard process (complete list of options).
    pythia.readString("SecondHard:generate = on");
    pythia.readString("SecondHard:TwoJets = on");
    //pythia.readString("SecondHard:PhotonAndJet = on");
    //pythia.readString("SecondHard:TwoPhotons = on");
    //pythia.readString("SecondHard:SingleGmZ = on");
    //pythia.readString("SecondHard:SingleW = on");
    //pythia.readString("SecondHard:TwoBJets = on");

    // Kinematics cuts, common for the two.
    pythia.readString("PhaseSpace:mHatMin = 40.");
    pythia.readString("PhaseSpace:pTHatMin = 20.");
  }

  // No parton or hadron level when creating LHEF. Limit output.
  if (outputLHEF) {
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:FSR = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");
  }
  pythia.readString("Next:numberShowLHA = 1");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");
  pythia.readString("Next:numberCount = 10000");
  pythia.readString("Check:nErrList = 2");

  // Optionally simplify run when comparing jets, for debug.
  if (findJets) {
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");
  }

  // Initialize.
  pythia.init();

  // Store initialization info in the LHAup object and write it out.
  if (outputLHEF) {
    myLHA.setInit();
    myLHA.initLHEF();
  }

  // Anti-kT jet finder. Exclude neutrinos from analysis.
  double etaMax   = 4.;
  double radius   = 0.7;
  double pTjetMin = 10.;
  int    nSel     = 2;
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  // Histogram.
  Hist wtViol( "violation weight distribution", 100, 0., 5.);
  Hist pTfirst("pT first collision",    100, 0., 400.);
  Hist pTsecond("pT second collision",  100, 0., 200.);
  Hist pTdiff("pT first-second collision", 100, -100., 300.);
  Hist nMult("number of multiparton interactions", 100, -0.5, 99.5);
  Hist bMore("b enhancement factor",    100, 0., 10.);
  Hist nChg("charged multiplicity", 100, -0.5, 999.5);
  Hist nFin("final multiplicity", 100, -0.5, 99.5);
  Hist nJet("jet multiplicity", 20, -0.5, 19.5);
  Hist pTj1("pT jet 1", 100, 0., 100.);
  Hist pTj2("pT jet 2", 100, 0., 100.);
  Hist pTj3("pT jet 3", 100, 0., 100.);
  Hist pTj4("pT jet 4", 100, 0., 100.);

  // Generate events.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    pythia.next();

    // Check if at end of the file.
    if (pythia.info.atEndOfFile()) break;

    // Check scales.
    if (iEvent < 10) cout << " process scales = " << fixed << setprecision(3) 
         << pythia.process.scale() << " and " << pythia.process.scaleSecond()
         << " and event scales = " << event.scale() << " and " 
         << event.scaleSecond() << endl;

    // Store event info in the LHAup object and write it out.
    if (outputLHEF) {
      myLHA.setEvent();
      myLHA.eventLHEF();
    }

    // Histogram event weight above unity.
    double wt = pythia.info.weight();
    if (wt > 1.) wtViol.fill( wt );

    // Histogram pT.
    double pT1 = pythia.info.pTMPI(0);
    double pT2 = pythia.info.pTMPI(1);
    pTfirst.fill( pT1 );
    pTsecond.fill( pT2 );
    pTdiff.fill( pT1 - pT2 );

    // Histogram multiparton interactions
    double nMPI = pythia.info.nMPI();
    nMult.fill( nMPI );
    bMore.fill( pythia.info.enhanceMPI() );

    // Histogram charged and final multiplicity.
    int nCharged = 0;
    int nFinal   = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      if (event[i].isCharged()) ++nCharged;
      ++nFinal;
    }
    nChg.fill( nCharged );
    nFin.fill( nFinal );

    // Analyze anti-kT jet properties and store result.
    if (findJets) {
      slowJet.analyze( pythia.event );
      int nJ = slowJet.sizeJet();
      nJet.fill( nJ );
      if (nJ > 0) pTj1.fill( slowJet.pT(0) );
      if (nJ > 1) pTj2.fill( slowJet.pT(1) );
      if (nJ > 2) pTj3.fill( slowJet.pT(2) );
      if (nJ > 3) pTj4.fill( slowJet.pT(3) );
    }

  }

  // Compare full statistics listing with what is set in info.
  pythia.stat();
  cout << scientific << setprecision(3) << "\n From pythia.info: sigma = "
       << pythia.info.sigmaGen() << " +- " << pythia.info.sigmaErr()
       << endl;

  // Update the cross section info and write endtag.
  if (outputLHEF) {
    myLHA.updateSigma();
    myLHA.closeLHEF(true);
  }

  // Print histograms.
  cout << wtViol << pTfirst << pTsecond << pTdiff << nMult << bMore << nChg 
       << nFin;
  if (findJets) cout << nJet << pTj1 << pTj2 << pTj3 << pTj4;

  // Done.
  return 0;
}
