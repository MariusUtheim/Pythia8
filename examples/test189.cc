// test189.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Particle production dependence on colour flow topology with
// different shower options.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 100000;

  // Process: 0 = only q q' -> q q', 1 = minbias, 2 = hard QCD.
  // For pure q q' -> q q' follow instructions in src/SigmaQCD.cc,
  // method double Sigma2qq2qq::sigmaHat().
  int iProc = 0;
  // pTmin for hard QCD processes.
  double pTmin = 45.;
  // pTmax, mMin and large scattering angles for q q' -> q q'.
  double pTmax = 50.;
  double mMin  = 500.;
  bool largeTheta = true;
  // Select large scattering angles
  // Include MPI or not.
  bool useMPI = false;

  // Histogram borders.
  double nchMax    = 199.;
  double deltaYmax = 8.;

  // Global histogram.
  Hist nchDeltaYRathist("<nch(Delta y)> old/new", 40, 0., deltaYmax);

  // Loop over shower options: 0 = none, 1 = old, 2 = new dipole.
  for (int modeShower = 0; modeShower < 3; ++modeShower) {

    // Generator. Shorthand for event.
    Pythia pythia;
    Event& event = pythia.event;

    // Set up incoming beams for 13 TeV LHC.
    pythia.readString("Beams:eCM = 13000.");

    // Set up QCD processes.
    if (iProc == 0) pythia.readString("HardQCD:qq2qq = on");
    if (iProc == 1) pythia.readString("SoftQCD:nonDiffractive = on");
    if (iProc == 2) pythia.readString("HardQCD:all = on");

    // Phase space cuts to probe forward-backward asymmetries.
    if (iProc == 0) pythia.settings.parm("PhaseSpace:mHatMin", mMin);
    if (iProc != 1) pythia.settings.parm("PhaseSpace:pTHatMin", pTmin);
    if (iProc == 0) pythia.settings.parm("PhaseSpace:pTHatMax", pTmax);
    if (iProc == 0 && largeTheta) pythia.settings.parm("PhaseSpace:Q2Min",
      0.5 * mMin * mMin);

    // MPI or not.
    if (!useMPI) pythia.readString("PartonLevel:MPI = off");

    // Reduce output.
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberCount = 100000");

    // Choose among shower options.
    if (modeShower == 0) {
      pythia.readString("PartonLevel:ISR = off");
      pythia.readString("PartonLevel:FSR = off");
    } else if (modeShower == 2) {
      pythia.readString("SpaceShower:dipoleRecoil = on");
    }

    // Initialize.
    pythia.init();

    // Histograms.
    Hist thetahist("theta*", 100, 0., M_PI);
    Hist deltaYhist("Delta y between partons", 40, 0., deltaYmax);
    Hist nchDeltaYhist("<nch(Delta y)>",       40, 0., deltaYmax);
    Hist nchhist("n_charged, inclusive",                100, -1., nchMax);
    Hist nchlohist("n_charged, low colour scattering",  100, -1., nchMax);
    Hist nchhihist("n_charged, high colour scattering", 100, -1., nchMax);

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Find charged multiplicity.
      int nch = 0;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && event[i].isCharged()) ++nch;
      nchhist.fill( nch );

      // Histogram multiplicity by rapidity separation between jets.
      double deltaY = abs( event[5].y() - event[6].y() );
      deltaYhist.fill( deltaY );
      nchDeltaYhist.fill( deltaY, nch );

      // Identify correct colour topologies.
      if ( event[3].id() > 0 && event[3].id() < 6 && event[4].id() > 0
        && event[4].id() < 6 && event[4].id() != event[3].id() ) {
        int iCol3 = 4;
        if (event[3].id() > 0 && event[5].col() == event[3].col()) iCol3 = 5;
        if (event[3].id() < 0 && event[5].acol() == event[3].acol()) iCol3 = 5;
        if (event[3].id() > 0 && event[6].col() == event[3].col()) iCol3 = 6;
        if (event[3].id() < 0 && event[6].acol() == event[3].acol()) iCol3 = 6;
        if (iCol3 != 4) {

          // Find hard-process scattering angle. Histogram multiplicity.
          RotBstMatrix Mrb;
          Mrb.toCMframe( event[3].p(), event[4].p() );
          Vec4 pOut = event[iCol3].p();
          pOut.rotbst( Mrb );
          double theta = pOut.theta();
          thetahist.fill( theta );
          if (theta < 0.5 * M_PI) nchlohist.fill( nch );
          else                    nchhihist.fill( nch );

        // End of colour topologies studies and event loop.
        }
      }
    }

    // Statistics and histograms.
    pythia.stat();
    nchDeltaYhist /= deltaYhist;
    cout << nchhist << deltaYhist << nchDeltaYhist
         << thetahist << nchlohist << nchhihist;

    // End loop over modeShower. Comparisons between shower options.
    if (modeShower == 1) nchDeltaYRathist  = nchDeltaYhist;
    if (modeShower == 2) nchDeltaYRathist /= nchDeltaYhist;
  }
  cout << nchDeltaYRathist;

  // Done.
  return 0;
}
