// test121.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Christine O. Rasmussen.

// The y, pT, x_Pomeron and t distributions for forward Z bosons at the LHC,
// within the hard diffraction framework for an inclusive event sample.
// Tests the impact of successive requirements.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Type of run (jets or Z) and number of events.
  bool isJets  = true;
  int nEvent   = 100000;

  // Create Pythia instance. Shorthand for event and info.
  Pythia pythia;
  Event& event = pythia.event;
  Info&  info  = pythia.info;

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
  pythia.readString("Diffraction:doHard = on");
  pythia.readString("Diffraction:sampleType = 1");
  pythia.readString("Diffraction:PomFlux = 5");
  pythia.readString("PDF:PomSet = 6");

  // Simplify printout. Print some diffractive events.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:showScaleAndVertex = off");
  int nPrintDiff = 2;
  int iPrintDiff = 0;

  // Switch off hadronization, since not used here.
  pythia.readString("HadronLevel:all = off");

  // Initialize.
  pythia.init();

  // Collect information on the number of diffractive events
  int nDiffA        = 0;
  int nDiffB        = 0;
  int nReducedDiffA = 0;
  int nReducedDiffB = 0;

  // Histograms.
  Hist y0("dN/dy inclusive",                 100, -5.,   5.);
  Hist y1("dN/dy after PDF selection",       100, -5.,   5.);
  Hist y2("dN/dy after MPI selection",       100, -5.,   5.);
  Hist yD1("dN/dy_diff after PDF selection", 100, -5.,   5.);
  Hist yD2("dN/dy_diff after MPI selection", 100, -5.,   5.);
  Hist pT0("dN/dpT inclusive",               100,  0., 100.);
  Hist pT1("dN/dpT after PDF selection",     100,  0., 100.);
  Hist pT2("dN/dpT after MPI selection",     100,  0., 100.);
  Hist pT3("dN/dpT after reweighting",       100,  0., 100.);
  Hist pT4("dN/dpT after/before reweighting",100,  0., 100.);
  Hist xP1("dN/dxPom after PDF selection",   100,  0.,   1.);
  Hist xP2("dN/dxPom after MPI selection",   100,  0.,   1.);
  Hist tP1("dN/dt after PDF selection",      100, -2.,   0.);
  Hist tP2("dN/dt after MPI selection",      100, -2.,   0.);
  Hist xI1("dN/dx_i/P after PDF selection",  100,  0.,   1.);
  Hist xI2("dN/dx_i/P after MPI selection",  100,  0.,   1.);
  Hist bI0("b_impact inclusive",             100,  0.,   5.);
  Hist bI1("b_impact after PDF selection",   100,  0.,   5.);
  Hist bI2("b_impact after MPI selection",   100,  0.,   5.);
  Hist xA0("dN/dx inclusive",                100,  0.,   1.);
  Hist xD1("dN/dx_i/P after PDF sel",        100,  0.,   1.);
  Hist xR1("dN/dx_rec after PDF sel",        100,  0.,   1.);
  Hist xD2("dN/dx_i/P after MPI sel",        100,  0.,   1.);
  Hist xR2("dN/dx_rec after MPI sel",        100,  0.,   1.);
  Hist xpA0("<x>(pT) inclusive",             100,  0., 100.);
  Hist xpD1("<x_i/P>(pT) after PDF sel",     100,  0., 100.);
  Hist xpR1("<x_rec>(pT) after PDF sel",     100,  0., 100.);
  Hist xpD2("<x_i/P>(pT) after MPI sel",     100,  0., 100.);
  Hist xpR2("<x_rec>(pT) after MPI sel",     100,  0., 100.);

  // Begin event loop. Generate event; skip if generation failed.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Hard process properties for jets or Z0.
    double yH  = info.y();
    double pTH = info.pTHat();
    double bI  = info.bMPI();
    if (!isJets) {
      int iZ = 0;
      for (int i = 0; i < event.size(); ++i) if (event[i].id() == 23) iZ = i;
      pTH = event[iZ].pT();
    }

    // Histogram inclusive properties.
    y0.fill( yH );
    pT0.fill( pTH );
    bI0.fill( bI );
    xA0.fill( info.x1pdf(), 0.5 );
    xA0.fill( info.x2pdf(), 0.5 );
    xpA0.fill( pTH, 0.5 * (info.x1pdf() + info.x2pdf()) );

    // Find diffractive events before MPI check.
    if (info.isHardDiffractive()) {
      bool isHDA = info.isHardDiffractiveA();

      // Properties of diffractive events.
      double yHD  = (isHDA) ? yH : -yH;
      double xPom = (isHDA) ? info.xPomeronB() : info.xPomeronA();
      double tPom = (isHDA) ? info.tPomeronB() : info.tPomeronA();
      double xPDF = (isHDA) ? info.x2pdf() : info.x1pdf();
      double xInP = xPDF/xPom;
      double xRec = (isHDA) ? info.x1pdf() : info.x2pdf();

      // Event weight from reweighted pdf.
      double wtX  = 1. / xInP;

      // Histogram properties before MPI check.
      if (isHDA) ++nDiffA; else ++nDiffB;
      y1.fill( yH );
      pT1.fill( pTH );
      bI1.fill( bI );
      yD1.fill( yHD );
      xP1.fill( xPom );
      tP1.fill( tPom );
      xI1.fill( xInP );
      xD1.fill( xInP );
      xR1.fill( xRec );
      xpD1.fill( pTH, xInP );
      xpR1.fill( pTH, xRec );

      // Find diffractive events passing MPI check.
      if (info.nMPI() == 1) {

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

        // Histogram properties after MPI check.
        if ( isHDA) ++nReducedDiffA; else ++nReducedDiffB;
        y2.fill( yH );
        pT2.fill( pTH );
        pT3.fill( pTH, wtX );
        bI2.fill( bI );
        yD2.fill( yHD );
        xP2.fill( xPom );
        tP2.fill( tPom );
        xI2.fill( xInP );
        xD2.fill( xInP );
        xR2.fill( xRec );
        xpD2.fill( pTH, xInP );
        xpR2.fill( pTH, xRec );
      }
    }

  // End of event loop. Statistics on event generation.
  }
  pythia.stat();

  // Statistics on diffraction.
  cout << "\n Side A is MPI-unchecked diffractive : " << nDiffA << endl;
  cout << " Side A is MPI-checked diffractive   : " << nReducedDiffA << endl;
  cout << " Side B is MPI-unchecked diffractive : " << nDiffB << endl;
  cout << " Side B is MPI-checked diffractive   : " << nReducedDiffB << endl;
  cout << " Total MPI-unchecked diffractive events : " << fixed
       << setprecision(2) << (nDiffA + nDiffB) / double(nEvent) * 100.
       << "%" << endl;
  cout << " Total MPI-checked diffractive events : "
       << (nReducedDiffA + nReducedDiffB) / double(nEvent) * 100.
       << "%" << endl;

  // Histograms.
  pT4 = pT3 / pT2;
  xpA0 /= pT0;
  xpD1 /= pT1;
  xpR1 /= pT1;
  xpD2 /= pT2;
  xpR2 /= pT2;
  cout << y0 << y1 << y2 << yD1 << yD2 << pT0 << pT1 << pT2 << pT3 << pT4
       << xP1 << xP2 << xI1 << xI2 << bI0 << bI1 << bI2
       << xA0 << xD1 << xR1 << xD2 << xR2
       << xpA0 << xpD1 << xpR1 << xpD2 << xpR2;
  //cout << tP1 << tP2;

  return 0;
}
