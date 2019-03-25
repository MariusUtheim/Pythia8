// test232.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug of dispaced vertices with beam vertex spread, for Felipe Rojas.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator. Shorthand for event record.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in  data from test232.cmnd.
  pythia.readFile( "test232.cmnd");

  // Extract data to be used in main program. Set counters.
  int nAbort  = pythia.mode("Main:timesAllowErrors");
  int iAbort  = 0;
  int idChi0   = 1000022;
  int idChiP   = 1000024;
  int nPrint   = 3;

  // Initialize generator.
  pythia.init();

  // Book histogram.
  Hist vPrim("vertex for primary chi0", 100, 1e-5, 1e5, true);
  Hist vSeco("vertex for secondary chi0", 100, 1e-5, 1e5, true);

  // Begin infinite event loop - to be exited at end of file.
  for (int iEvent = 0; ; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) {

      // Leave event loop if at end of file.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    // Print first few events.
    if (iEvent < nPrint) {
      pythia.LHAeventList();
      pythia.process.list(true);
      pythia.event.list(true);
    }

    // Search for final chi0.
    for (int i = 1; i < event.size(); ++i)
    if (event[i].id() == idChi0 && event[i].isFinal()) {
      Vec4 vProd = event[i].vProd();
      double dProd = max( 1.01e-5, vProd.pAbs() );
      int idMot = event[event[i].mother1()].idAbs();
      if      (idMot == idChi0) vPrim.fill( dProd);
      else if (idMot == idChiP) vSeco.fill( dProd);
      else cout << " Warning: unexpected chi0 mother " << idMot << endl;
    }

  // End of event loop.
  }

  // Give statistics. Print histograms.
  pythia.stat();
  cout << vPrim << vSeco;
  HistPlot hpl("plot232");
  hpl.frame( "fig232", "Displacement of chi0 production vertex",
    "distance (mm)");
  hpl.add( vPrim, "h,red", "primary chi0");
  hpl.add( vSeco, "h,blue", "chi0 from chi+- decay");
  hpl.plot();

  // Done.
  return 0;
}
