// test107.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// pT spectrum in diffeent multiplicity bins.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events. Pseudorapidity cut.
  int nEvent = 1000;
  double etaMax = 2.5;

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Histograms.
  Hist mult("charged multiplicity", 100, -0.5, 499.5);
  Hist pTall("dn_ch/dpT all", 100, 0., 10.);
  Hist pT1("(1/n_ch) dn_ch/dpT n_ch < 20", 100, 0., 10.);
  Hist pT2("(1/n_ch) dn_ch/dpT 40 < n_ch < 60", 100, 0., 10.);
  Hist pT3("(1/n_ch) dn_ch/dpT n_ch < 90", 100, 0., 10.);
  Hist pTnow("dn_ch/dpT current event", 100, 0., 10.);
  int n1 = 0, n2 = 0, n3 = 0;

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of final charged particles and their pT spectrum.
    int nCh = 0;
    pTnow.null();
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged()
       && abs(pythia.event[i].eta()) < etaMax) {
        ++nCh;
        pTnow.fill( pythia.event[i].pT() );
    }

    // Fill histograms for multiplicity and normalized pT distribution.
    mult.fill( nCh );
    pTall += pTnow;
    if (nCh < 20) {++n1; pT1 += pTnow / nCh;}
    if (nCh > 40 && nCh < 600) {++n2; pT2 += pTnow / nCh;}
    if (nCh > 90) {++n3; pT3 += pTnow / nCh;}

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  pTall *= 10. / nEvent;
  pT1 *= 10. / n1;
  pT2 *= 10. / n2;
  pT3 *= 10. / n3;
  cout << mult << pTall << pT1 << pT2 << pT3;

  // Done.
  return 0;
}
