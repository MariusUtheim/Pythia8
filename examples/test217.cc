// test217.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study the interference contribution to elastic scattering.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. t range. Histogram array.
  int    nEvent  = 1000000;
  int    nAbort  = 5;
  double tAbsMin = 2e-5;
  double tAbsMid = 1e-3;
  double tAbsMax = 0.2;
  Hist tEl1[12], tEl2[12], tRat[12];

  // Loop over cases. Book histograms. Generator.
  for (int ica = 0; ica < 12; ++ica) {
    tEl1[ica].book( "  ", 100, tAbsMin, tAbsMax, true);
    tEl2[ica].book( "  ", 100, tAbsMin, tAbsMax, true);
    tRat[ica].book( "  ", 100, tAbsMin, tAbsMax, true);
    Pythia pythia;
    cout << "\n Now running case number " << ica << endl;

    // Common settings.
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("SoftQCD:elastic = on");
    pythia.settings.parm("SigmaElastic:tAbsMin", tAbsMin);
    pythia.readString("SigmaElastic:rho = 0.15");
    pythia.readString("Next:numberCount = 0");

    // Case-by-case settings. 0-3: SaS; 4-7: ABMST; 8-11: RPP.
    // 0-1, 4-5, 8-9 pp, rest ppbar. Even without, odd with Coulomb.
    if      (ica < 4)   pythia.readString("SigmaTotal:mode = 1");
    else if (ica < 8)   pythia.readString("SigmaTotal:mode = 3");
    else                pythia.readString("SigmaTotal:mode = 4");
    if ((ica/2)%2 == 1) pythia.readString("Beams:idB = -2212");
    if (ica%2 == 1)     pythia.readString("SigmaElastic:Coulomb = on");

    // Initialize.
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

      // Study t distribution of elastic events.
      double tAbs = abs(pythia.info.tHat());
      tEl1[ica].fill(tAbs);
      if (tAbs > tAbsMid) tEl2[ica].fill(tAbs);

    // End of event loop. Normalize histograms. End case loop.
    }
    pythia.stat();
    double sig = pythia.info.sigmaGen();
    tEl1[ica] *= sig / nEvent;
    tEl2[ica] *= sig / nEvent;
  }

  // Ratios.
  tRat[0] = tEl1[1]/tEl1[3];
  tRat[1] = tEl1[5]/tEl1[7];
  tRat[2] = tEl1[9]/tEl1[11];
  for (int j = 0; j < 6; ++j)
    tRat[3 + j] = (tEl2[2 * j + 1] - tEl2[2 * j]) / tEl2[2 * j];

  // Histograms.
  HistPlot hpl("test217plot");
  hpl.frame( "out217plot", "Elastic $|t|$ spectrum N+C", "$|t|$ (GeV$^2$)",
    "$\\mathrm{d}N / \\mathrm{d}\\,\\log |t|$");
  hpl.add( tEl1[1], "-,red", "SaS pp");
  hpl.add( tEl1[3], "--,red", "SaS ppbar");
  hpl.add( tEl1[5], "-,green", "ABMST pp");
  hpl.add( tEl1[7], "--,green", "ABMST ppbar");
  hpl.add( tEl1[9], "-,blue", "RPP pp");
  hpl.add( tEl1[11], "--,blue", "RPP ppbar");
  hpl.plot();
  hpl.frame( "", "Elastic $|t|$ ratio pp/ppbar", "$|t|$ (GeV$^2$)",
    "$\\mathrm{d}N / \\mathrm{d}|t|$ ratio");
  hpl.add( tRat[0] , "-,red", "SaS pp/ppbar");
  hpl.add( tRat[1] , "-,green", "ABMST pp/ppbar");
  hpl.add( tRat[2] , "-,blue", "RPP pp/ppbar");
  hpl.plot();
  hpl.frame( "", "Elastic $|t|$ ratio (N+C)/N $-$ 1", "$|t|$ (GeV$^2$)",
    "$\\mathrm{d}N / \\mathrm{d}|t|$ ratio");
  hpl.add( tRat[3], "-,red", "SaS pp" );
  hpl.add( tRat[4], "--,red", "SaS ppbar" );
  hpl.add( tRat[5], "-,green", "ABMST pp" );
  hpl.add( tRat[6], "--,green", "ABMST ppbar" );
  hpl.add( tRat[7], "-,blue", "RPP pp" );
  hpl.add( tRat[8], "--,blue", "Rpp ppbar" );
  hpl.plot();

  // Done.
  return 0;
}
