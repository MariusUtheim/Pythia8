// test107.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Helical event structure.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate.
  int nEvent = 10000;

  // Histograms.
  Hist Seta( "S_eta(xi) with phi",      50, 0., 5.);
  Hist Setax("S_eta(xi) phi = 0 or pi", 50, 0., 5.);
  Hist Seta0("S_eta(xi) phi = 0",       50, 0., 5.);
  Hist SetaS("S_eta(xi) swapped eta",   50, 0., 5.);
  Hist SetaRx("S_eta(xi) phi = 0 or pi/ phi = 0", 50, 0., 5.);
  Hist SetaR0("S_eta(xi) with phi / phi = 0", 50, 0., 5.);
  Hist SetaRR("S_eta(xi) phi = 0 or pi/ with phi", 50, 0., 5.);

  // Variables for analysis and statistics.
  int nAcc = 0;
  vector<double> etaJ, etaS, phiJ, phiX;
  double xiVal[50];
  for (int k = 0; k < 50; ++k) xiVal[k] = 0.1 * (k + 0.5);
  complex<double> phaseVal[50], phaseValx[50], phaseVal0[50], phaseValS[50];
  double rfac, phaseValH;

  // Generator; shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Select process.
  pythia.readString("SoftQCD:nonDiffractive = on");

  // Restrict output.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize.
  pythia.init();

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Skip if failure.
    if (!pythia.next()) continue;

    // Reset variables.
    bool isAcc = true;
    etaJ.clear();
    etaS.clear();
    phiJ.clear();
    phiX.clear();
    for (int k = 0; k < 50; ++k) {
      phaseVal[k]  = 0.;
      phaseValx[k] = 0.;
      phaseVal0[k] = 0.;
      phaseValS[k] = 0.;
    }

    // Pick up acceptable tracks.
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].isCharged()
//      && abs(event[i].eta()) < 2.5 && event[i].pT() > 0.1) {
//      && abs(event[i].eta()) < 4.0 && event[i].pT() > 0.1) {
      && abs(event[i].eta()) < 1.5 && event[i].pT() > 0.1) {
      if (event[i].pT() > 10.) {isAcc = false; break;}
      etaJ.push_back(event[i].eta());
      etaS.push_back(event[i].eta());
      phiJ.push_back(event[i].phi());
      phiX.push_back( (event[i].px() > 0.) ? 0. : M_PI );
    }

    // Skip ahead for unacceptable events.
    int mult = etaJ.size();
    if (mult < 6) isAcc = false;
    if (!isAcc) continue;

    // Rearranged rapidities.
    for (int i = 0; i < mult; ++i) {
      int j = min( int( mult * pythia.rndm.flat() ), mult - 1);
      swap( etaS[i], etaS[j] );
    }

    // Sum up phases for different xi values.
    ++nAcc;
    for (int i = 0; i < mult; ++i)
    for (int k = 0; k < 50; ++k) {
      phaseVal[k]  += exp( complex<double>(0., xiVal[k] * etaJ[i] - phiJ[i]) );
      phaseValx[k] += exp( complex<double>(0., xiVal[k] * etaJ[i] - phiX[i]) );
      phaseVal0[k] += exp( complex<double>(0., xiVal[k] * etaJ[i]) );
      phaseValS[k] += exp( complex<double>(0., xiVal[k] * etaS[i] - phiJ[i]) );
    }

    // Histogram values for S_eta(xi) - 1.
    rfac = 1. / mult;
    for (int k = 0; k < 50; ++k) {
      phaseValH = rfac * pow2(abs(phaseVal[k]));
      Seta.fill( xiVal[k], phaseValH - 1.);
      phaseValH = rfac * pow2(abs(phaseValx[k]));
      Setax.fill( xiVal[k], phaseValH - 1.);
      phaseValH = rfac * pow2(abs(phaseVal0[k]));
      Seta0.fill( xiVal[k], phaseValH - 1.);
      phaseValH = rfac * pow2(abs(phaseValS[k]));
      SetaS.fill( xiVal[k], phaseValH - 1.);
    }

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Print histograms.
  Seta /= nAcc;
  Setax /= nAcc;
  Seta0 /= nAcc;
  SetaS /= nAcc;
  SetaRx = Setax / Seta0;
  SetaR0 = Seta / Seta0;
  SetaRR = Setax / Seta;
  cout << Seta << Setax << Seta0 << SetaS << SetaRx << SetaR0 << SetaRR;

  // Done.
  return 0;
}
