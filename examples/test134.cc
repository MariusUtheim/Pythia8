// test134.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Christine O. Rasmussen.

// Various distributions for jet or Z boson production at the LHC,
// within the hard diffraction framework for an inclusive event sample.
// Tests the impact of successive requirements.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Type of run (jets or Z) and number of events without/with MPI check.
  bool isJets   = true;
  int nEvent1   = 4000;
  int nEvent2   = 10 * nEvent1;

  // PomFlux: 1 = SaS, 2 = Bruni-Ingelman, 3 = Berger-Streng, 4 = DL,
  //          5 = MBR, 6 = H1 Fit A, 7 = H1 Fit B.
  int pomFlux   = 7;

  // PomSet: 1 = simple, 2 = pi0, 3 = H1 Fit A NLO, 4 = H1 Fit B NLO,
  //         5 = H1 Jets NLO, 6 = H1 Fit B LO.
  int pomSet    = 6;

  // bSelHard: 1 = same b as pp, 2 = sqrt(old b), 3 = new b.
  int bSelHard  = 2;

  // Histograms.
  Hist y0("dN/dy inclusive",                 100, -5.,   5.);
  Hist y1("dN/dy after PDF selection",       100, -5.,   5.);
  Hist y2("dN/dy after MPI selection",       100, -5.,   5.);
  Hist yD1("dN/dy_diff after PDF selection", 100, -5.,   5.);
  Hist yD2("dN/dy_diff after MPI selection", 100, -5.,   5.);
  Hist pT0("dN/dpT inclusive",               100,  0., 100.);
  Hist pT1("dN/dpT after PDF selection",     100,  0., 100.);
  Hist pT2("dN/dpT after MPI selection",     100,  0., 100.);
  Hist bI0("b_impact inclusive",             100,  0.,   5.);
  Hist bI1("b_impact after PDF selection",   100,  0.,   5.);
  Hist bI2("b_impact after MPI selection",   100,  0.,   5.);
  Hist bI3("b_impact for Pom-p after MPI",   100,  0.,   5.);
  Hist bI4("P_{nMPI=1}(b_impact) inclusive", 100,  0.,   5.);
  Hist xI1("dN/dx_i/P after PDF selection",  100,  0.,   1.);
  Hist xI2("dN/dx_i/P after MPI selection",  100,  0.,   1.);
  Hist xP1("dN/dxPom after PDF selection",   100,  0.,   1.);
  Hist xP2("dN/dxPom after MPI selection",   100,  0.,   1.);
  Hist xA0("dN/dx inclusive (both sides)",   100,  0.,   1.);
  Hist xT1("dN/d(x_i/P * xPom) after PDF",   100,  0.,   1.);
  Hist xT2("dN/d(x_i/P * xPom) after MPI",   100,  0.,   1.);
  Hist xR1("dN/dx_rec after PDF sel",        100,  0.,   1.);
  Hist xR2("dN/dx_rec after MPI sel",        100,  0.,   1.);
  Hist tP1("dN/dt after PDF selection",      100, -2.,   0.);
  Hist tP2("dN/dt after MPI selection",      100, -2.,   0.);

  // Loop over diffraction without/with MPI check.
  for (int iMPI = 0; iMPI < 3; ++iMPI) {
    int nEvent = (iMPI != 2) ? nEvent1 : nEvent2;

    // Create Pythia instance. Shorthand for event and info.
    Pythia pythia;
    Event& event   = pythia.event;
    Info&  info    = pythia.info;

    // Set it up to generate jets or Z's at 8 TeV.
    pythia.readString("Beams:eCM = 8000.");
    if (isJets) {
      pythia.readString("HardQCD:all = on");
      pythia.readString("PhaseSpace:pTHatMin = 20.");
    } else {
      pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
      pythia.readString("23:mMin = 70.");
      pythia.readString("23:mMax = 110.");
    }

    // Setup of diffractive framework.
    if (iMPI == 0) pythia.readString("Diffraction:doHard = off");
    else           pythia.readString("Diffraction:doHard = on");
    pythia.settings.mode("Diffraction:sampleType", iMPI);
    pythia.settings.mode("Diffraction:PomFlux", pomFlux);
    pythia.settings.mode("PDF:PomSet", pomSet);
    pythia.settings.mode("Diffraction:bSelHard", bSelHard);

    // Simplify printout. Print some diffractive events.
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:showScaleAndVertex = off");
    pythia.readString("Next:numberCount = 0");
    int nPrintDiff = 0;
    int iPrintDiff = 0;

    // Switch off hadronization, since not used here.
    pythia.readString("HadronLevel:all = off");

    // Initialize.
    pythia.init();

    // Collect information on the number of diffractive events.
    int nAcc   = 0;
    int nDiffA = 0;
    int nDiffB = 0;

    // Begin event loop. Generate event; skip if generation failed.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Hard process properties for jets or Z0.
      double yH  = info.y();
      double pTH = info.pTHat();
      double bI  = (iMPI < 2) ? info.bMPI() : info.bMPIold();
      if (!isJets) {
        int iZ = 0;
        for (int i = 0; i < event.size(); ++i)
          if (event[i].id() == 23) iZ = i;
        pTH = event[iZ].pT();
      }

      // Histogram inclusive properties.
      if (iMPI == 0) {
        ++nAcc;
        y0.fill( yH );
        pT0.fill( pTH );
        bI0.fill( bI );
        if (info.nMPI() == 1) bI4.fill( bI );
        xA0.fill( info.x1pdf(), 0.5 );
        xA0.fill( info.x2pdf(), 0.5 );
      }

      // Find diffractive events.
      if (info.isHardDiffractive()) {
        bool isHDA  = info.isHardDiffractiveA();
        if (isHDA) ++nDiffA;
        else       ++nDiffB;

        // Print info on first few diffractive events passing MPI check.
        if (iPrintDiff++ < nPrintDiff) {
          info.list();
          pythia.process.list();
          cout << " isHardDiff = " << info.isHardDiffractiveA() << "  "
               << info.isHardDiffractiveB() << endl << scientific;
          cout << " xPom = " << info.xPomeronA() << "  "
               << info.xPomeronB() << endl;
          cout << " xPDF = " << info.x1pdf() << "  "
               << info.x2pdf() << endl;
        }

        // Properties of diffractive events.
        double yHD  = (isHDA) ? yH : -yH;
        double xPom = (isHDA) ? info.xPomeronB() : info.xPomeronA();
        double xPDF = (isHDA) ? info.x2pdf() : info.x1pdf();
        double xInP = xPDF/xPom;
        double xRec = (isHDA) ? info.x1pdf() : info.x2pdf();
        double tPom = (isHDA) ? info.tPomeronB() : info.tPomeronA();

        // Histogram properties without MPI check.
        if (iMPI == 1) {
          ++nAcc;
          y1.fill( yH );
          yD1.fill( yHD );
          pT1.fill( pTH );
          bI1.fill( bI );
          xI1.fill( xInP );
          xP1.fill( xPom );
          xT1.fill( xInP * xPom );
          xR1.fill( xRec );
          tP1.fill( tPom );

        // Histogram properties after MPI check.
        } else {
          ++nAcc;
          y2.fill( yH );
          yD2.fill( yHD );
          pT2.fill( pTH );
          bI2.fill( bI );
          bI3.fill( info.bMPI() );
          xI2.fill( xInP );
          xP2.fill( xPom );
          xT2.fill( xInP * xPom );
          xR2.fill( xRec );
          tP2.fill( tPom );
        }

      // End of hard diffraction part and of event loop.
      }
    }

    // Statistics on event generation.
    pythia.stat();

    // Statistics on diffraction.
    double topc = 100. / double(nEvent);
    cout << "\n Side A is diffractive : " << nDiffA << endl;
    cout << " Side B is diffractive : " << nDiffB << endl;
    cout << " Total diffractive events : " << fixed
       << setprecision(2) << topc * (nDiffA + nDiffB) << "%" << endl;

    // Rescale histograms.
    if (iMPI == 0) {
      y0   *= 10.  / nAcc;
      pT0  *= 1.   / nAcc;
      bI4  /= bI0;
      bI0  *= 20.  / nAcc;
      xA0  *= 100. / nAcc;
    } else if (iMPI == 1) {
      y1   *= 10.  / nAcc;
      yD1  *= 10.  / nAcc;
      pT1  *= 1.   / nAcc;
      bI1  *= 20.  / nAcc;
      xI1  *= 100. / nAcc;
      xP1  *= 100. / nAcc;
      xT1  *= 100. / nAcc;
      xR1  *= 100. / nAcc;
      tP1  *= 50.  / nAcc;
    } else {
      y2   *= 10.  / nAcc;
      yD2  *= 10.  / nAcc;
      pT2  *= 1.   / nAcc;
      bI2  *= 20.  / nAcc;
      bI3  *= 20.  / nAcc;
      xI2  *= 100. / nAcc;
      xP2  *= 100. / nAcc;
      xT2  *= 100. / nAcc;
      xR2  *= 100. / nAcc;
      tP2  *= 50.  / nAcc;
    }

  // End loop over diffraction without/with MPI check.
  }

  // Histograms.
  cout << y0 << y1 << y2 << yD1 << yD2 << pT0 << pT1 << pT2
       << bI0 << bI1 << bI2 << bI3 << bI4 << xI1 << xI2 << xP1 << xP2
       << xA0 << xT1 << xT2 << xR1 << xR2 << tP1 << tP2;

  return 0;
}
