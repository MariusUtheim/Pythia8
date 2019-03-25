// test218.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study the hadronic pT spectrum by multiplicity class.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Include diffraction and colour reconnection or not.
  int nEvent  = 100000;
  int nAbort  = 5;
  bool doDiff = true;
  bool doCR   = true;

  // Book histograms.
  int nm[8] = {0};
  Hist nchH( "charged multiplicity", 100, -0.5, 99.5);
  Hist nchdH( "charged multiplicity diffractive", 100, -0.5, 99.5);
  Hist pTnH[8];
  for (int im = 0; im < 8; ++im) pTnH[im].book( "  ", 100, 0.1, 10., true);

  // Generator setup.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");
  if (doDiff) pythia.readString("SoftQCD:inelastic = on");
  else        pythia.readString("SoftQCD:nondiffractive = on");
  if (!doCR)  pythia.readString("ColourReconnection:reconnect = off");
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

    // Special classification diffractive.
    bool isDiff = (pythia.info.code() != 101);
    ++nm[0];
    if (isDiff) ++nm[1];

    // Hadronic charged multiplicity in |eta| < 0.8 and its pT spectrum.
    int nch = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()
      && event[i].isHadron() && abs(event[i].eta()) < 0.8) {
      ++nch;
      pTnH[0].fill( event[i].pT() );
      if (isDiff) pTnH[1].fill( event[i].pT() );
    }
    nchH.fill( nch );
    if (isDiff) nchdH.fill( nch );

    // Classify events by multiplicity above.
    int im = 2;
    if (nch >  5) im = 3;
    if (nch > 10) im = 4;
    if (nch > 20) im = 5;
    if (nch > 30) im = 6;
    if (nch > 40) im = 7;
    ++nm[im];

    // Fill pT spectra by multiplicity class.
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()
      && event[i].isHadron() && abs(event[i].eta()) < 0.8) {
      pTnH[im].fill( event[i].pT() );
    }

  // End of event loop. Normalize histograms. End case loop.
  }
  pythia.stat();
  nchH /= nm[0];
  nchdH /= nm[0];
  for (int im = 0; im < 8; ++im) pTnH[im] /= nm[im];
  for (int im = 2; im < 8; ++im) pTnH[im] /= pTnH[0];

  // Histograms.
  HistPlot hpl("test218plot");
  hpl.frame( "out218plot", "Hadronic charged multiplicity in $|\\eta| < 0.8$",
    "$n_h^{\\pm}$", "$\\mathrm{d}P / \\mathrm{d}n_h^{\\pm}$");
  hpl.add( nchH, ",red", "all");
  if (doDiff) hpl.add( nchdH, ",blue", "diffractive");
  hpl.plot();
  hpl.frame( "", "Hadronic $p_{\\perp}$ spectrum in $|\\eta| < 0.8$",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}P / \\mathrm{d}\\,\\log(p_{\\perp})$");
  hpl.add( pTnH[0] , ",red", "all");
  if (doDiff) hpl.add( pTnH[1] , ",blue", "diffractive");
  hpl.plot();
  hpl.frame( "", "$p_{\\perp}$ ratios by multiplicity class",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}P / \\mathrm{d}p_{\\perp}$ ratio");
  hpl.add( pTnH[2], "", "$n_h^{\\pm} \\leq 5$");
  hpl.add( pTnH[3], "", "$6 \\leq n_h^{\\pm} \\leq 10$");
  hpl.add( pTnH[4], "", "$11 \\leq n_h^{\\pm} \\leq 20$");
  hpl.add( pTnH[5], "", "$21 \\leq n_h^{\\pm} \\leq 30$");
  hpl.add( pTnH[6], "", "$31 \\leq n_h^{\\pm} \\leq 40$");
  hpl.add( pTnH[7], "", "$n_h^{\\pm} \\geq 41$");
  hpl.plot(true);

  // Done.
  return 0;
}
