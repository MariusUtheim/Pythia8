// test905.cc
// ISR pT of gamma*/Z0 for varying cm energy or mass.
// Exercise 5.2.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events, energies, mass ranges.
  int nEvent = 10000;
  double eCM[4] = {2000., 4000., 8000., 16000.};
  double mMin[4] = {  90., 22.5, 45., 180.};
  double mMax[4] = { 110., 27.5, 55., 220.};

  // Loop over CM energies and Z mass range.
  for (int iCase = 0; iCase < 7; ++iCase) {
    int iEcm = (iCase < 4) ? iCase : 2;
    int iMgZ = (iCase < 4) ? 0 : iCase - 3;

    // Set up Pythia; code common for cases.
    Pythia pythia;
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("ProcessLevel:resonanceDecays = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");

    // Set up cm energy and gamma*/Z0 mass range.
    pythia.settings.parm("Beams:eCM", eCM[iEcm]);
    pythia.particleData.mMin(23, mMin[iMgZ]);
    pythia.particleData.mMax(23, mMax[iMgZ]);

    // Initialize Pythia.
    pythia.init();
    Event& event = pythia.event;

    // Histogram and statistics for Z pT spectrum.
    Hist ZpT("pT of Z0", 100, 0.,100.);
    double avgPT = 0.;

    // Generate events.
    for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if ( !pythia.next() ) continue;

      // Find and histogram final Z0.
      int iZ = 0;
      for (int i = 0; i < event.size(); ++i) if (event[i].id() == 23) iZ = i;
      ZpT.fill( event[iZ].pT() );
      avgPT +=  event[iZ].pT();

    // End of event and case loops. Print histogram and statistics.
    }
    cout << ZpT;
    cout << fixed << setprecision(3) << "\n Average Z pT is "
         << avgPT / nEvent << " for Ecm = " << eCM[iEcm]
         << " and mass around " << 0.5 * (mMin[iMgZ] + mMax[iMgZ]) << endl;
  }

  return 0;
}
