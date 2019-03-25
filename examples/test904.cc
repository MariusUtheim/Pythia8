// test904.cc
// Parton emission multiplicities at a few e+e- energies.
// Exercise 3.5

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events, energies, statistics.
  int nEvent = 10000;
  double eCM[4] = {25., 50., 100., 200.};
  int nFin[4] = {0};
  int nEmt[4] = {0};

  // Loop over CM energies.
  for (int iEcm = 0; iEcm < 4; ++iEcm) {

    // Set up Pythia.
    Pythia pythia;
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.settings.parm("Beams:eCM", eCM[iEcm]);
    pythia.readString("PDF:lepton = off");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("HadronLevel:all = off");
    pythia.readString("23:onMode = off");
    //pythia.readString("23:onIfAny = 1 2 3 4 5");
    pythia.readString("23:onIfAny = 1 2");
    pythia.init();
    Event& event = pythia.event;

    // Generate events.
    for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if ( !pythia.next() ) continue;

      // Loop through particles: final and gluons emitted.
      for (int i = 0; i < event.size(); ++i) {
        if (event[i].isFinal() ) ++nFin[iEcm];
        if (event[i].id() == 21 && event[event[i].mother1()].idAbs() < 6)
          ++nEmt[iEcm];

      // End of particle, event and energy loops.
      }
    }
  }

  // Print statistics at the four energies.
  for (int iEcm = 0; iEcm < 4; ++iEcm) {
    cout << fixed << setprecision(3) << " At Ecm = " << eCM[iEcm]
         << " n_final is " <<  double(nFin[iEcm])/double(nEvent)
         << " and n_q/qbar->g is " << double(nEmt[iEcm])/double(nEvent)
         << endl;
  }

  return 0;
}
