// test910.cc
// Size of edge pseudorapidity gap.
// Exercise 9.3.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events. QCD or top.
  int nEvent = 10000;

  // Histograms for pseudorapidity gap.
  Hist etaGapND( "nondiffractive", 100, 0., 10.);
  Hist etaGapSD( "single diffractive", 100, 0., 10.);

  // Set up and initialize Pythia for ND or SD.
  for (int iDiff = 0; iDiff < 2; ++iDiff) {
    Pythia pythia;
    pythia.readString("Beams:Ecm = 13000.");
    if (iDiff == 0) pythia.readString("SoftQCD:nonDiffractive = on");
    else            pythia.readString("SoftQCD:singleDiffractive = on");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.init();
    Event& event = pythia.event;

    // Generate events.
    for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if ( !pythia.next() ) continue;

      // Find largest and smallest eta inside detector.
      int nInside = 0;
      double etaMax = -5.;
      double etaMin = 5.;
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].pT() > 0.2) {
        double etaNow = event[i].eta();
        if (abs(etaNow) < 5.) {
          ++nInside;
          if (etaNow > etaMax) etaMax = etaNow;
          if (etaNow < etaMin) etaMin = etaNow;
        }
      }

      // Fill histogram with largest edge pseudorapidity gap.
      if (nInside > 0) {
        double etaGapNow = max( 5. - etaMax, 5. + etaMin);
        if (iDiff == 0) etaGapND.fill( etaGapNow );
        else            etaGapSD.fill( etaGapNow );
      }

    // End of event loop. End of ND/SD loop.
    }
  }

  // Print histograms.
  etaGapND *= 10./nEvent;
  etaGapSD *= 10./nEvent;
  cout << etaGapND << etaGapSD;
  HistPlot hpl("test910plot");
  hpl.frame( "out910plot", "Edge rapidity gap size",
    "$\\Delta y_{\\mathrm{gap}}$", "Probability");
  hpl.add( etaGapND );
  hpl.add( etaGapSD );
  hpl.plot( true);

  // Done.
  return 0;
}
