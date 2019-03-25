// test137.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Number of partons above threshold.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Key Switches.
  int nEvent = 10000;
  bool doHard = true;

  // Set up soft or hard QCD.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 2760.");
  if (doHard) {
    pythia.readString("PartonLevel:all = on");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 10.");
    pythia.readString("Tune:pp = 5");
  } else {
    pythia.readString("SoftQCD:nonDiffractive = on");
  }
  pythia.readString("HadronLevel:all = off");
  pythia.init();

  // Histogram and statistics.
  Hist pTqH( "pT spectrum for quarks", 100, 0., 100.);
  Hist pTgH( "pT spectrum for gluons", 100, 0., 100.);
  int nq10 = 0;
  int nq20 = 0;
  int ng10 = 0;
  int ng20 = 0;
  double pTq[100];
  double pTg[100];
  for(int i = 0; i < 100; i++) {
    pTq[i] = 0.0;
    pTg[i] = 0.0;
  }
  double pTnow;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find final-state partons. Fill statistics.
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      pTnow = event[i].pT();
      if (event[i].idAbs() < 6) {
        pTqH.fill( pTnow );
        if (pTnow > 10.) ++nq10;
        if (pTnow > 20.) ++nq20;
      } else if (event[i].id() == 21) {
        pTgH.fill( pTnow );
        if (pTnow > 10.) ++ng10;
        if (pTnow > 20.) ++ng20;
      }
      if (pTnow < 100.0) {
        int pTidx = (int)(pTnow);
        if (event[i].idAbs() < 6) pTq[pTidx]++;
        else if (event[i].id() == 21)  pTg[pTidx]++;
      }
    }
  }

  // Statistics.
  pythia.stat();
  cout << pTqH;
  cout << "\n number of quarks above 10 GeV = " << nq10
       << " and above 20 GeV = " << nq20 << endl;
  cout << pTgH;
  cout << "\n number of gluons above 10 GeV = " << ng10
       << " and above 20 GeV =  " << ng20 << endl;
  ofstream quark_dis("quark_spectrum.dat");
  ofstream gluon_dis("gluon_spectrum.dat");
  for(int i = 0; i < 100; i++) {
    quark_dis << scientific << setw(18) << setprecision(8)
              << (i+0.5) << "   " << (double)pTq[i]/(double)nEvent << endl;
    gluon_dis << scientific << setw(18) << setprecision(8)
              << (i+0.5) << "   " << (double)pTg[i]/(double)nEvent << endl;
  }
  quark_dis.close();
  gluon_dis.close();

  return 0;
}
