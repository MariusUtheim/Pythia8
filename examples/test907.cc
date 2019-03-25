// test907.cc
// Nondiffractive "min bias" event properties.
// Exercise 6.3 and 6.4.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 10000;

  // Set up and initialize Pythia.
  Pythia pythia;
  pythia.readString("Beams:Ecm = 8000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  //pythia.readString("PartonLevel:MPI = off");
  pythia.readString("ColourReconnection:reconnect = off");
  pythia.init();
  Event& event = pythia.event;

  // Histogram and statistics on particle production.
  Hist nchD("charged multiplicity", 100, -1., 399.);
  Hist etaD("charged peudorapidity", 100, -10., 10.);
  Hist pTD("charged pT spectrum", 100, 0., 10.);
  Hist pTnchD("charged <pT>(n)", 100, -1., 399.);
  double yMin[9] = {0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.};
  double yMax[9] = {1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.};
  double nFs[9]  = {0.};
  double nBs[9]  = {0.};
  double nF2s[9] = {0.};
  double nB2s[9] = {0.};
  double nFBs[9] = {0.};

  // Generate events.
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if ( !pythia.next() ) continue;

    // Reset statistics for current event.
    int nch     = 0;
    double pTch = 0.;
    int nF[9]   = {0};
    int nB[9]   = {0};

    // Inclusive charged particle properties in the event.
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].isCharged()) {
       ++nch;
       etaD.fill( event[i].eta() );
       pTD.fill( event[i].pT() );
       pTch += event[i].pT();

       // Particle multiplicity in rapidity bins.
       double yNow = event[i].y();
       if (yNow > 0.) {for (int j = 0; j < 9; ++j)
         if (yNow > yMin[j] && yNow < yMax[j]) ++nF[j];}
       else           {for (int j = 0; j < 9; ++j)
         if (yNow > -yMax[j] && yNow < -yMin[j]) ++nB[j];}
    }

    // Fill/save integrated properites of events.
    nchD.fill( nch );
    pTnchD.fill( nch, pTch / nch );
    for (int j = 0; j < 9; ++j) {
      nFs[j]  += nF[j];
      nBs[j]  += nB[j];
      nF2s[j] += nF[j] * nF[j];
      nB2s[j] += nB[j] * nB[j];
      nFBs[j] += nF[j] * nB[j];
    }
  }

  // Normalize and print histograms and forward-backward correlations.
  pTnchD /= nchD;
  nchD *= 0.5 / nEvent;
  etaD *=  5. / nEvent;
  pTD  *= 10. / nEvent;
  cout << nchD << etaD << pTD << pTnchD << endl;
  for (int j = 0; j < 9; ++j) {
    double nFavg2  = pow2( (nFs[j] + nBs[j]) / (2. * nEvent) );
    double nF2avg = (nF2s[j] + nB2s[j]) / (2. * nEvent);
    double nFBavg = nFBs[j] / nEvent;
    double rhoFB  = (nFBavg - nFavg2) / (nF2avg - nFavg2);
    cout << " rho_FB(Delta-y = " << fixed << setprecision(2) << double(j)
         << ") = " << setprecision(3) << rhoFB << endl;
  }

  return 0;
}
