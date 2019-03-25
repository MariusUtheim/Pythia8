// test216.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study the interference contribution to elastic scattering.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. t range.
  int    nEvent  = 1000000;
  int    nAbort  = 5;
  double tAbsMin = 2e-5;
  double tAbsMax = 0.2;

  // Book histograms.
  Hist tElLog0("SaS pp",      100, tAbsMin, tAbsMax, true);
  Hist tElLog1("SaS ppbar",   100, tAbsMin, tAbsMax, true);
  Hist tElLog2("ABMST pp",    100, tAbsMin, tAbsMax, true);
  Hist tElLog3("ABMST ppbar", 100, tAbsMin, tAbsMax, true);
  Hist tElLog4("RPP pp",      100, tAbsMin, tAbsMax, true);
  Hist tElLog5("RPP ppbar",   100, tAbsMin, tAbsMax, true);
  Hist tElLog6("SaS pp N+C",  100, tAbsMin, tAbsMax, true);
  Hist tElLog7("SaS pp N",    100, tAbsMin, tAbsMax, true);
  Hist tElLog8("SaS ppbar N+C", 100, tAbsMin, tAbsMax, true);
  Hist tElLog9("SaS ppbar N",   100, tAbsMin, tAbsMax, true);

  // Loop over cases. Generator.
  for (int ica = 0; ica < 10; ++ica) {
    Pythia pythia;

    // Settings to be used in the main program.
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("SoftQCD:elastic = on");
    if (ica != 7 && ica != 9) pythia.readString("SigmaElastic:Coulomb = on");
    pythia.settings.parm("SigmaElastic:tAbsMin", tAbsMin);
    if (ica > 5) pythia.settings.parm("SigmaElastic:tAbsMin", 1e-3);
    if (ica == 2 || ica == 3) pythia.readString("SigmaTotal:mode = 3");
    if (ica == 4 || ica == 5) pythia.readString("SigmaTotal:mode = 4");
    if (ica == 1 || ica == 3 || ica == 5 || ica > 7)
      pythia.readString("Beams:idB = -2212");
    pythia.readString("SigmaElastic:rho = 0.15");
    pythia.readString("Next:numberCount = 0");

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
      if (ica == 0) tElLog0.fill(tAbs);
      if (ica == 1) tElLog1.fill(tAbs);
      if (ica == 2) tElLog2.fill(tAbs);
      if (ica == 3) tElLog3.fill(tAbs);
      if (ica == 4) tElLog4.fill(tAbs);
      if (ica == 5) tElLog5.fill(tAbs);
      if (ica == 6 && tAbs > 1e-3) tElLog6.fill(tAbs);
      if (ica == 7 && tAbs > 1e-3) tElLog7.fill(tAbs);
      if (ica == 8 && tAbs > 1e-3) tElLog8.fill(tAbs);
      if (ica == 9 && tAbs > 1e-3) tElLog9.fill(tAbs);

    // End of event loop. Normalize histogram. End case loop.
    }
    pythia.stat();
    double sig = pythia.info.sigmaGen();
    if (ica == 0) tElLog0 *= sig / nEvent;
    if (ica == 1) tElLog1 *= sig / nEvent;
    if (ica == 2) tElLog2 *= sig / nEvent;
    if (ica == 3) tElLog3 *= sig / nEvent;
    if (ica == 4) tElLog4 *= sig / nEvent;
    if (ica == 5) tElLog5 *= sig / nEvent;
    if (ica == 6) tElLog6 *= sig / nEvent;
    if (ica == 7) tElLog7 *= sig / nEvent;
    if (ica == 8) tElLog8 *= sig / nEvent;
    if (ica == 9) tElLog9 *= sig / nEvent;
  }

  // Ratios.
  Hist tElRat0("SaS   pp/ppbar", 100, tAbsMin, tAbsMax, true);
  Hist tElRat1("ABMST pp/ppbar", 100, tAbsMin, tAbsMax, true);
  Hist tElRat2("RPP   pp/ppbar", 100, tAbsMin, tAbsMax, true);
  Hist tElRat3("SaS pp (N+C)/N - 1", 100, tAbsMin, tAbsMax, true);
  Hist tElRat4("SaS ppbar (N+C)/N - 1", 100, tAbsMin, tAbsMax, true);
  tElRat0 = tElLog0 / tElLog1;
  tElRat1 = tElLog2 / tElLog3;
  tElRat2 = tElLog4 / tElLog5;
  tElRat3 = (tElLog6 - tElLog7) / tElLog7;
  tElRat4 = (tElLog8 - tElLog9) / tElLog9;

  // Histograms.
  cout << tElLog0 << tElLog1;
  HistPlot hpl("test216plot");
  hpl.frame( "out216plot", "Elastic $|t|$ spectrum", "$|t|$ (GeV$^2$)",
    "$(1/|t|) \\mathrm{d}N / \\mathrm{d}|t|$");
  hpl.add( tElLog0, "-");
  hpl.add( tElLog1, "-");
  hpl.add( tElLog2, "-");
  hpl.add( tElLog3, "-");
  hpl.add( tElLog4, "-");
  hpl.add( tElLog5, "-");
  hpl.plot();
  hpl.frame( "", "Elastic $|t|$ spectrum", "$|t|$ (GeV$^2$)",
    "$(1/|t|) \\mathrm{d}N / \\mathrm{d}|t|$");
  hpl.add( tElRat0, "-");
  hpl.add( tElRat1, "-");
  hpl.add( tElRat2, "-");
  hpl.plot();
  hpl.frame( "", "Elastic $|t|$ spectrum (N+C)/N - 1", "$|t|$ (GeV$^2$)",
    "$(1/|t|) \\mathrm{d}N / \\mathrm{d}|t|$");
  hpl.add( tElRat3, "-");
  hpl.add( tElRat4, "-");
  hpl.plot();

  // Done.
  return 0;
}
