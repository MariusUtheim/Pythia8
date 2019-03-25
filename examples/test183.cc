// test183.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Jet cross section as a function of pTmin.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Common parameters of integration.
  int npT = 50;
  double pTmin = 0.8;
  double pTmax = 20.;
  int nEvt = 100000;

  // Run 1 = SppS, 2 = Tevatron, 3 = LHC.
  for (int run = 1; run < 4; ++run) {
  double eCM[4] = { 0., 630., 1960., 13000.};

  // Output file.
  string fileOut[4] = {" ", "sigmasps.data", "sigmatev.data", "sigmalhc.data"};
  ofstream strmOut(fileOut[run].c_str());

  // Get total cross section.
  Pythia pythiaTot;
  if (run < 3) pythiaTot.readString("Beams:idB = -2212");
  pythiaTot.settings.parm("Beams:eCM", eCM[run]);
  pythiaTot.readString("SoftQCD:all = on");
  pythiaTot.readString("PartonLevel:all = off");
  pythiaTot.readString("Print:quiet = on");
  pythiaTot.init();
  for (int iEvt = 0; iEvt < nEvt; ++iEvt) pythiaTot.next();
  double sigmaTot = pythiaTot.info.sigmaGen();

  // Loop over pTmin values.
  Pythia pythia;
  for (int ipT = 0; ipT < npT; ++ ipT) {
    double pTnow = pTmin * pow( pTmax / pTmin, ipT / (npT - 1.) );

    // Generator. Process selection.
    if (run < 3) pythia.readString("Beams:idB = -2212");
    pythia.settings.parm("Beams:eCM", eCM[run]);
    pythia.readString("HardQCD:all = on");
    pythia.settings.forceParm("PhaseSpace:pThatMin", pTnow);
    pythia.settings.forceParm("PhaseSpace:pTHatMinDiverge", pTmin);

    // Switch off unnecessary parts. Initialize
    pythia.readString("PartonLevel:all = off");
    pythia.readString("Print:quiet = on");
    pythia.init();

    // Event loop.
    for (int iEvt = 0; iEvt < nEvt; ++iEvt) pythia.next();

    // Statistics.
    double sigma = pythia.info.sigmaGen();
    cout << scientific << setprecision(4) << setw(12) << pTnow
         << setw(12) << sigma << setw(12) << sigmaTot << endl;
    strmOut << scientific << setprecision(4) << setw(12) << pTnow
         << setw(12) << sigma << setw(12) << sigmaTot << endl;
  }

  // End loop over runs.
  }

  return 0;
}
