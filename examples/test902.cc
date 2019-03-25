// test902.cc
// Pseudorapidity spectra.
// Exercise 2.7

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 10000;

  // Set up Pythia.
  Pythia pythia;
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();
  Event& event = pythia.event;

  // Histograms.
  Hist ypi("y_pi" , 100, -10.,10.);
  Hist yK("y_K" ,   100, -10.,10.);
  Hist yp("y_p" ,   100, -10.,10.);
  Hist etapi("eta_pi" , 100, -10.,10.);
  Hist etaK("eta_K" ,   100, -10.,10.);
  Hist etap("eta_p" ,   100, -10.,10.);

  // Generate events.
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if ( !pythia.next() ) continue;

    // Loop through final particles.
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() ) {
      if (event[i].idAbs() == 211) {
        ypi.fill( event[i].y() );
        etapi.fill( event[i].eta() );
      } else if (event[i].idAbs() == 321) {
        yK.fill( event[i].y() );
        etaK.fill( event[i].eta() );
      } else if (event[i].idAbs() == 2212) {
        yp.fill( event[i].y() );
        etap.fill( event[i].eta() );
      }
    }
  }

  // Print histograms.
  cout << ypi << yK << yp << etapi << etaK << etap;

  return 0;
}
