// test222.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study flavour composition as a function of multiplicity.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 10000000;
  int nAbort  = 5;
  int nTrig   = 0;

  // Choice of model: 0 = default, 1 = Jesper, 2 = Nadine.
  int model = 2;

  // Hadron species to study: 0 = ch, 1 = pi+-, 2 = K0S, 3 = p,
  // 4 = Lambda, 5 = Xi, 6 = Omega.
  int idSel[7] = { 0, 211, 310, 2212, 3122, 3312, 3334};
  long nSum[7] = { 0, 0, 0, 0, 0, 0, 0};
  int nNow[7];

  // Book histograms.
  Hist nHad[7];
  for (int im = 0; im < 7; ++im) nHad[im].book( "  ", 24, 0.5, 24.5);

  // Common generator setup.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 7000.");
  pythia.readString("SoftQCD:inelastic = on");

  // Jesper and Peter's CR model.
  if (model == 1) {
    pythia.readString("ColourReconnection:mode = 1");
    pythia.readString("ColourReconnection:allowDoubleJunRem = off");
    pythia.readString("BeamRemnants:remnantMode = 1");
    pythia.readString("BeamRemnants:saturation = 5");
    pythia.readString("StringPT:sigma = 0.335");
    pythia.readString("StringZ:aLund = 0.36");
    pythia.readString("StringZ:bLund = 0.56");
    pythia.readString("StringFlav:probQQtoQ = 0.078");
    pythia.readString("StringFlav:ProbStoUD = 0.217");
    pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");

  // Nadine's model.
  } else if (model == 2) {
    pythia.readString("StringPT:thermalModel = on");
    pythia.readString("StringPT:temperature = 0.21");
    pythia.readString("StringPT:closePacking = on");
    pythia.readString("ColourReconnection:range = 1.1");
    pythia.readString("MultipartonInteractions:pT0Ref = 2.5");
    pythia.readString("HadronLevel:HadronScatter = on");
    //pythia.readString("");
  }

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

    // Events should have at least one charged particle in |eta| < 1.
    int nch1 = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()
      && event[i].isCharged() && abs(event[i].eta()) < 1.) ++nch1;
    if (nch1 == 0) continue;
    ++nTrig;

    // Charged and other multiplicity in |eta| < 0.5.
    for (int im = 0; im < 7; ++im) nNow[im] = 0;
    for (int i = 0; i < event.size(); ++i) if (abs(event[i].eta()) < 0.5) {
      if (event[i].isFinal() && event[i].isCharged()) ++nNow[0];
      int idAbs = event[i].idAbs();
      for (int im = 1; im < 7; ++im) if (idAbs == idSel[im]) ++nNow[im];
    }
    nHad[0].fill( nNow[0], 1. );
    for (int im = 1; im < 7; ++im) nHad[im].fill( nNow[0], nNow[im] );

    // End of event loop. Normalize histograms.
    for (int im = 0; im < 7; ++im) nSum[im] += nNow[im];
  }
  pythia.stat();
  nHad[0] /= nTrig;
  for (int im = 1; im < 7; ++im) nHad[im] /= double(nTrig) * double(nSum[im]);
  for (int im = 2; im < 7; ++im) nHad[im] /= nHad[1];

  // Plot histograms.
  HistPlot hpl("test223plot");
  hpl.frame( "out223plot", "Charged multiplicity in $|\\eta| < 0.5$",
    "$n_{\\mathrm{ch}}$", "$\\mathrm{d}P / \\mathrm{d}n_{\\mathrm{h}}$");
  hpl.add( nHad[0], "", "charged");
  hpl.plot();
  hpl.frame( "", "(h/$\\pi$) normalized",
    "$n_{\\mathrm{ch}}$ in $|\\eta| < 0.5$",
    "(h/$\\pi$)($n_{\\mathrm{ch}}$) / (ditto all $n_{\\mathrm{ch}}$)");
  hpl.add( nHad[2], "-,green", "$K^0_S$");
  hpl.add( nHad[3], "-,red", "p");
  hpl.add( nHad[4], "-,blue", "$\\Lambda$");
  hpl.add( nHad[5], "-,magenta", "$\\Xi$");
  hpl.add( nHad[6], "-,cyan", "$\\Omega$");
  hpl.plot();
  hpl.frame( "out223plotX", "(h/$\\pi$) normalized",
    "$n_{\\mathrm{ch}}$ in $|\\eta| < 0.5$",
    "(h/$\\pi$)($n_{\\mathrm{ch}}$) / (ditto all $n_{\\mathrm{ch}}$)");
  hpl.add( nHad[2], "-,green", "$K^0_S$");
  hpl.add( nHad[3], "-,red", "p");
  hpl.add( nHad[4], "-,blue", "$\\Lambda$");
  hpl.add( nHad[5], "-,magenta", "$\\Xi$");
  hpl.plot();

  // Done.
  return 0;
}
