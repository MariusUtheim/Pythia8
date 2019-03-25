// test207.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Compare dipole recoil schemes for colour singlet exchange process.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent    = 100000;

  // Histograms.
  Hist yJet12old("y for jet 1 and 2, old", 100, -5., 5.);
  Hist yJet3old("y for jet 3, old", 100, -5., 5.);
  Hist yJet12new("y for jet 1 and 2, new", 100, -5., 5.);
  Hist yJet3new("y for jet 3, new", 100, -5., 5.);
  Hist yJet12rat("y for jet 1 and 2, new/old", 100, -5., 5.);
  Hist yJet3rat("y for jet 3, new/old", 100, -5., 5.);
  Hist yShiftold("y(jet - parton), old", 100, -2., 2.);
  Hist yShiftnew("y(jet - parton), new", 100, -2., 2.);

  // Set up anti-kT jet finder. (SlowJet is frontend to FJcore.)
  double etaMax   = 4.5;
  double radius   = 0.4;
  double pTjetMin = 20.;
  int    nSel     = 2;
  int    massSet  = 1;
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);

  // Loop over old and new recoil schemes.
  for (int iRec = 0; iRec < 2; ++iRec) {
    Hist& yJet12 = (iRec == 0) ? yJet12old : yJet12new;
    Hist& yJet3  = (iRec == 0) ? yJet3old  : yJet3new;
    Hist& yShift = (iRec == 0) ? yShiftold : yShiftnew;

    // Generator. Shorthand for event.
    Pythia pythia;
    Event& event = pythia.event;

    // Set up WW-fusion process q q -> q q H, with Higgs stable.
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("HiggsSM:ff2Hff(t:WW) = on");
    pythia.readString("25:mayDecay = off");

    // Switch on new ISR recoil scheme.
    if (iRec == 1) pythia.readString("SpaceShower:dipoleRecoil = on");

    // Simplify output. Initialize.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.init();

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Forcibly remove Higgs from analysis.
      for (int i = 0; i < event.size(); ++i)
        if (event[i].id() == 25) event[i].statusNeg();

      // Analyze jet properties.
      slowJet.analyze( event );
      int nJet = slowJet.sizeJet();

      // Fill jet rapidities.
      if (nJet > 0) yJet12.fill( slowJet.y(0) );
      if (nJet > 1) yJet12.fill( slowJet.y(1) );
      if (nJet > 2) yJet3.fill( slowJet.y(2) );

      // Correlate jets with scattered quarks.
      Vec4 pQ1 = event[6].p();
      Vec4 pQ2 = event[7].p();
      for (int iJet = 0; iJet < min( 2, nJet); ++iJet) {
        Vec4 pJet = slowJet.p(iJet);
        double dR1 = RRapPhi( pJet, pQ1);
        double dR2 = RRapPhi( pJet, pQ2);
        if ( min( dR1, dR2) < 2.) {
	  int iQ = (dR1 < dR2) ? 6 : 7;
	  double dy = slowJet.y(iJet) - event[iQ].y();
	  if (event[iQ].y() < 0.) dy = -dy;
	  yShift.fill( dy );
        }
      }

    // End of event loop.
    }

    // Statistics. End of recoil scheme loop.
    pythia.stat();
  }

  // Histograms.
  yJet12rat = yJet12new / yJet12old;
  yJet3rat = yJet3new / yJet3old;
  cout << yJet12old << yJet12new << yJet12rat << yJet3old << yJet3new
       << yJet3rat << yShiftold << yShiftnew;

  // Done.
  return 0;
}
