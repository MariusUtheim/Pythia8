// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("PhaseSpace:mHatMin = 80.");
  pythia.readString("PhaseSpace:mHatMax = 120.");
  pythia.readString("23:mayDecay = off");
  pythia.init();
  Hist pTZ("dN/dpTZ", 100, 0., 100.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Loop over particles in event. Find last Z0 copy. Fill its pT.
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (event[i].id() == 23) iZ = i;
    pTZ.fill( event[iZ].pT() );

    // Decay the Z0.
    if (!event[iZ].isFinal()) cout << " Error: final Z not stable " << endl;
    double mZ    = event[iZ].m();
    double mJpsi = pythia.particleData.m0(443);
    double pAbs  = 0.5 * (mZ*mZ - mJpsi*mJpsi) / mZ;
    Vec4 pJpsi(  0., 0.,  pAbs, 0.5 * (mZ*mZ + mJpsi*mJpsi) / mZ);
    Vec4 pGamma( 0., 0., -pAbs, pAbs);
    double theta = acos( 2. * pythia.rndm.flat() - 1.);
    double phi   = 2. * M_PI * pythia.rndm.flat();
    pJpsi.rot(  theta, phi);
    pGamma.rot( theta, phi);
    pJpsi.bst(  event[iZ].p() );
    pGamma.bst( event[iZ].p() );
    int iDau1 = event.append( 443, 91, iZ, 0, 0, 0, 0, 0, pJpsi, mJpsi);
    int iDau2 = event.append(  22, 91, iZ, 0, 0, 0, 0, 0, pGamma, 0.);
    event[iZ].daughters( iDau1, iDau2);
    event[iZ].statusNeg();
    pythia.moreDecays();

    // Illustrate modified event record.
    if (iEvent == 0) event.list();

  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << pTZ;
  return 0;
}
