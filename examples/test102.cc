// test102.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Number of soft jets around a hard q or g jet in Z + jet,
// for Mihoko Nojiri and Bryan Webber, December 2014.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events. Main selection.
  int nEvent      = 5000;
  double pTHatMin = 450.;
  double pTJetMin = 500.;
  int nListJets   = 1;

  // Loop over q and g jet.
  for (int iQG = 0; iQG < 2; ++iQG) {

    // Generator. Output restriction.
    Pythia pythia;
    pythia.readString("Next:numberShowEvent = 0");

    // Process selection.
    pythia.readString("Beams:eCM = 13000.");
    if (iQG == 0) pythia.readString("WeakBosonAndParton:qg2gmZq = on");
    else          pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
    pythia.settings.parm("PhaseSpace:pTHatMin", pTHatMin);

    // Low gamma*/Z0 masses are cut away, and Z0 put stable.
    pythia.readString("23:mMin = 80.");
    pythia.readString("23:mayDecay = off");

    // Special selections for ISR off and no MPI/FSR interleaving.
    //pythia.readString("PartonLevel:ISR = off");
    //pythia.readString("PartonLevel:MPI = off");
    //pythia.readString("TimeShower:interleave = off");

    // Initialization.
    pythia.init();

    // Anti-kT jet finder. Exclude neutrinos from analysis. Jet separation.
    double jetRadius = 0.4;
    double jetpTmin  = 20.;
    SlowJet slowJet( -1, jetRadius, jetpTmin, 10., 2, 1);
    double jetSepMax = 0.8;

    // Histograms.
    Hist nJetsTot("number of jets in total", 50, -0.5, 49.5);
    Hist pTJet0("pT hardest jet", 100, 0., 1000.);
    Hist pTJetI("pT all non-hardest jets", 100, 0., 1000.);
    Hist pTHatAcc("pT of the hard interaction, accepted", 100, 0., 1000.);
    Hist pTHatRej("pT of the hard interaction, rejected", 100, 0., 1000.);
    Hist nRingAcc("number of jets around hardest, accepted", 10, -0.5, 9.5);
    Hist nRingRej("number of jets around hardest, rejected", 10, -0.5, 9.5);

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Remove final Z0 from analysis.
      for (int i = 0; i < pythia.event.size(); ++i)
        if (pythia.event[i].id() == 23) pythia.event[i].statusNeg();

      // Find jets in event. List jets in first few events.
      slowJet. analyze( pythia.event );
      if (iEvent < nListJets) slowJet.list();

      // Fill inclusive jet distributions and above pT threshold.
      int nJets = slowJet.sizeJet();
      nJetsTot.fill( nJets );
      double pT0 = slowJet.pT(0);
      pTJet0.fill( pT0 );
      for (int i = 1; i < slowJet.sizeJet(); ++i)
        pTJetI.fill( slowJet.pT(i) );
      bool accept = (pT0 > pTJetMin);
      if (accept) pTHatAcc.fill( pythia.info.pTHat() );
      else        pTHatRej.fill( pythia.info.pTHat() );

      // Count number of jets in ring around hardest.
      int nRing = 0;
      for (int i = 1; i < slowJet.sizeJet(); ++i) {
        double dY   = slowJet.y(i) - slowJet.y(0);
        double dPhi = abs( slowJet.phi(i) - slowJet.phi(0) );
        if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
        double dR = sqrt( pow2(dY) + pow2(dPhi) );
        if (dR < jetSepMax) ++nRing;
      }
      if (accept) nRingAcc.fill( nRing );
      else        nRingRej.fill( nRing );


    // End of event loop. Statistics. Histograms.
    }
    pythia.stat();
    cout << nJetsTot << pTJet0 << pTJetI << pTHatAcc << pTHatRej
         << nRingAcc << nRingRej;

  // End of loop over q and g jets.
  }

  // Done.
  return 0;
}
