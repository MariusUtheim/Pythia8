// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// The pT_W spectrum at the LHC, shower vs. analytic suppression.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent = 10000;

  // Book histograms.
  Hist pTWshower("shower", 100, 0., 50.);
  Hist pTWme("matrix elements", 100, 0., 50.);
  Hist pTWmeSudno("ME + Sudakov, no alphas prefac", 100, 0., 50.);
  Hist pTWmeSudyes("ME + Sudakov, with alphas prefac", 100, 0., 50.);
  Hist pTWmeSudmid("ME + Sudakov, intermediate alphas", 100, 0., 50.);

  // Loop over W and W + 1 jet options. Create alpha_strong object.
  for (int nJet = 0; nJet < 2; ++nJet) {
    AlphaStrong alphaS;
    alphaS.init(0.130);

    // Generator. Process selection.
    Pythia pythia;
    pythia.readString("Beams:eCM = 7000.");
    if (nJet == 0) {
      pythia.readString("WeakSingleBoson:ffbar2W = on");
    } else {
      pythia.readString("WeakBosonAndParton:qqbar2Wg = on");
      pythia.readString("WeakBosonAndParton:qg2Wq = on");
      pythia.readString("PhaseSpace:pTHatMin = 0.5");
      pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");
    }
    pythia.readString("24:mMin = 60.");
    pythia.readString("24:mMax = 100.");

    // Simplify generation and initialize.
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PartonLevel:FSR = off");
    if (nJet == 1) pythia.readString("PartonLevel:ISR = off");
    pythia.readString("HadronLevel:all = off");
    pythia.readString("BeamRemnants:primordialKT = off");
    pythia.init();

    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Loop over particles in event. Find last W+- copy. Fill its pT.
      int iW = 0;
      for (int i = 0; i < pythia.event.size(); ++i)
        if (pythia.event[i].idAbs() == 24) iW = i;
      double pTW = pythia.event[iW].pT();
      if (nJet == 0) pTWshower.fill( pTW );
      else if (pTW > 2.) pTWme.fill( pTW );

      // Multiply by Sudakov and alpha_s factor for W + 1 jet.
      if (nJet == 1) {
        double mW = pythia.event[iW].m();
        double logRatio = 2. * log( pTW / mW );
        double w = exp( - 2. * ((4./3.) * alphaS.alphaS(pTW*pTW) / (2. * M_PI))
          * ( 0.5 * pow2(logRatio) + 1.5 * logRatio ) );
        pTWmeSudno.fill( pTW, w);
        w *= alphaS.alphaS( pTW*pTW ) / pythia.info.alphaS();
        pTWmeSudyes.fill( pTW, w);
        w = exp( - 2. * ((4./3.) * alphaS.alphaS(pTW * mW) / (2. * M_PI))
          * ( 0.5 * pow2(logRatio) + 1.5 * logRatio ) );
        pTWmeSudmid.fill( pTW, w);
      }

    // End of event loop. Statistics.
    }
    pythia.stat();

    // Normalize histograms. End of jets loop.
    double facNorm = 1e9 * pythia.info.sigmaGen() / (0.5 * nEvent);
    if (nJet == 0) pTWshower *= facNorm;
    if (nJet == 1) { pTWme *= facNorm; pTWmeSudno *= facNorm;
      pTWmeSudyes *= facNorm; pTWmeSudmid *= facNorm; }
  }

  // Print histograms.
  cout << pTWshower << pTWme << pTWmeSudno << pTWmeSudyes << pTWmeSudmid;

  // Write Python code to print the histograms in the same plot(s).
  HistPlot hpl("test209");
  hpl.frame( "WpT1", "W $p_{\\perp}$ spectrum", "$p_{\\perp}$ (GeV)",
    "d$\\sigma$/d$p_{\\perp}$ (pb/GeV)");
  hpl.add( pTWshower);
  hpl.add( pTWme, "-");
  //hpl.add( pTWmeSudno, ".");
  hpl.add( pTWmeSudyes, "--");
  //hpl.add( pTWmeSudmid, "-.");
  hpl.plot(true);
  hpl.frame( "WpT2", "W $p_{\\perp}$ spectrum", "$p_{\\perp}$ (GeV)",
    "d$\\sigma$/d$p_{\\perp}$ (pb/GeV)");
  hpl.add( pTWshower);
  hpl.add( pTWmeSudno, ".", "ME + Sud, fix $\\alpha_s$");
  hpl.add( pTWmeSudyes, "--", "ME + Sud + $\\alpha_s$ reweight");
  hpl.add( pTWmeSudmid, "-.", "ME + Sud(intermediate $\\alpha_s$");
  hpl.plot();

  return 0;
}
