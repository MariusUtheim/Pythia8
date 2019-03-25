// test201.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test program of diffractive scattering in various frameworks.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Choose pp/ppbar, CM energy and number of events.
  double eCM         = 10000.;
  int    nEvent      = 1000;

  // Histogram.
  Hist epsZero( "epsilon = 0",    80, 1., 10000., true);
  Hist epsPos(  "epsilon = 0.1",  80, 1., 10000., true);
  Hist epsNeg(  "epsilon = -0.1", 80, 1., 10000., true);

  // Loop over three epsilon values. Generator.
  for (int iEps = 0; iEps < 3; ++iEps) {
  Pythia pythia;

  // Set up run: diffractive cross sections.
  //pythia.readString("SoftQCD:all = on");
  //pythia.readString("SoftQCD:singleDiffractive = on");
  pythia.readString("SoftQCD:doubleDiffractive = on");
  //pythia.readString("SoftQCD:centralDiffractive = on");
  pythia.readString("SigmaTotal:zeroAXB = off");
  pythia.readString("SigmaTotal:sigmaAXB2TeV = 2.");
  pythia.readString("SigmaDiffractive:maxAXB = 5.");
  if (iEps == 1) pythia.readString("SigmaDiffractive:SaSepsilon = 0.1");
  if (iEps == 2) pythia.readString("SigmaDiffractive:SaSepsilon = -0.1");

  // Final common setup. Initialization.
  pythia.settings.parm("Beams:eCM", eCM);
  pythia.readString("PartonLevel:all = off");
  int nAbort = 5;
  pythia.init();

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Diffractive mass scale.
    double mNow = 1.;
    for (int i = 3; i < pythia.process.size(); ++i)
    if (pythia.process[i].m() > mNow) mNow = pythia.process[i].m();
    if (iEps == 0) epsZero.fill( mNow );
    if (iEps == 1) epsPos.fill( mNow );
    if (iEps == 2) epsNeg.fill( mNow );

  // End of event loop.
  }

  // Final statistics. Normalize and print histograms.
  pythia.stat();
  }
  cout << epsZero << epsPos << epsNeg;

  // Done.
  return 0;
}
