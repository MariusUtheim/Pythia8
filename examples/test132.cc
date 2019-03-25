// test132.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check the pShift and daughterListRecursive methods.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Generator.
  Pythia pythia;
  Rndm& rndm = pythia.rndm;
  double px, py, pz, e, m1, m2, m3, m4;
  double pMax = 2.;

  // Check that momentum shift works.
  for (int i = 0; i < 100; ++i) {

    // New and old masses.
    m1 = rndm.flat();
    m2 = rndm.flat();
    m3 = rndm.flat();
    m4 = rndm.flat();

    // Pick three-momenta at random.
    px = pMax * (2. * rndm.flat() - 1.);
    py = pMax * (2. * rndm.flat() - 1.);
    pz = pMax * (2. * rndm.flat() - 1.);
    e  = sqrt(px*px + py*py + pz*pz + m1*m1);
    Vec4 p1( px, py, pz, e);
    px = pMax * (2. * rndm.flat() - 1.);
    py = pMax * (2. * rndm.flat() - 1.);
    pz = pMax * (2. * rndm.flat() - 1.);
    e  = sqrt(px*px + py*py + pz*pz + m2*m2);
    Vec4 p2( px, py, pz, e);

    cout << fixed << setprecision(3) << "\n Shift from " << m1 << " and "
         << m2 << " to " << m3 << " and " << m4 << " for mHat = "
         << (p1 + p2).mCalc() << endl;
    cout << " pIn1  = " << p1 << " pIn2  = " << p2;

    // Shuffle momenta.
    bool isOK = pShift( p1, p2, m3, m4);
    if (isOK) {
      cout << " pOut1 = " << p1 << " pOut2 = " << p2;
      double err = abs(p1.mCalc() - m3) + abs(p2.mCalc() - m4);
      if (err > 1e-6) cout << scientific << " Error is " << err << endl;
    } else cout << " reshuffling failed!" << endl;
  }

  // Set up top production and decay.
  pythia.readString("Top:all = on");
  pythia.readString("HadronLevel:all = off");
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 20; ++iEvent) {
    if (!pythia.next()) continue;

    // List daughterListRecursive.
    pythia.process.list();
    for (int i = 1; i < pythia.process.size(); ++i) {
      vector<int> iDau = pythia.process[i].daughterListRecursive();
      cout << " i = " << i << " has daughters";
      for (int j = 0; j < int(iDau.size()); ++j) cout << " " << iDau[j];
      cout << endl;
    }

  }

  return 0;
}
