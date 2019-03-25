// test917.cc
// Generate fragmentation according to the Artru-Mennesier model.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//-----------------------------------------------------------------

int main() {

  // Number of events; CM energy; number of vertices per event.
  int    nEvent = 100000;
  double eCM    = 10.;
  int    nVtx   = 50;
  int    nAcc, nOrd;
  double pPVtx[50], pMVtx[50], pPAcc[50], pMAcc[50], pPOrd[52], pMOrd[52];
  double m2;
  bool   doSave;

  // Set up Pythia for use as random number generator.
  Pythia pythia;
  pythia.init();

  // Book histograms.
  Hist nAccH("number of accepted vertices", 50, -0.5, 49.5);
  Hist mAccH("mass of accepted hadrons", 100, 0., 10.);

  // Loop over events.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate the vertices at random inside full phase space.
    for (int iVtx = 0; iVtx < nVtx; ++iVtx) {
      pPVtx[iVtx] = eCM * pythia.rndm.flat();
      pMVtx[iVtx] = eCM * pythia.rndm.flat();
    }

    // Skip those that are in the forward lightcone of another.
    nAcc = 0;
    for (int iVtx = 0; iVtx < nVtx; ++iVtx) {
      doSave = true;
      for (int jVtx = 0; jVtx < nVtx; ++jVtx) if (jVtx != iVtx
      && pPVtx[iVtx] > pPVtx[jVtx] && pMVtx[iVtx] > pMVtx[jVtx]) {
        doSave = false;
        break;
      }

      // Add to list of accepted vertices.
      if (doSave) {
        pPAcc[nAcc] = pPVtx[iVtx];
        pMAcc[nAcc] = pMVtx[iVtx];
        ++nAcc;
      }
    }

    // Order in terms of decreasing pPlus, with endpoints added.
    pPOrd[0] = eCM;
    pMOrd[0] = 0.;
    for (int i = 0; i < nAcc; ++i) {
      pPOrd[i + 1] = pPAcc[i];
      pMOrd[i + 1] = pMAcc[i];
      for (int j = i; j > 0; --j) {
        if (pPOrd[j] > pPOrd[j + 1]) break;
        swap( pPOrd[j], pPOrd[j + 1]);
        swap( pMOrd[j], pMOrd[j + 1]);
      }
    }
    pPOrd[nAcc + 1] = 0.;
    pMOrd[nAcc + 1] = eCM;
    nOrd = nAcc + 2;

    // Print result for first few events.
    if (iEvent < 5) {
      cout << endl << endl << fixed << setprecision(3);
      for (int i = 0; i < nOrd; ++i) cout << setw(3) << i
        << setw(8) << pPOrd[i] << setw(8) << pMOrd[i] << endl;
    }

    // Histogram number and mass.
    nAccH.fill( nAcc );
    for (int i = 1; i < nOrd; ++i) {
      m2 = (pPOrd[i - 1] - pPOrd[i]) * (pMOrd[i] - pMOrd[i - 1]);
      mAccH.fill( sqrt(m2) );
    }

  // End of event loop. Print Histograms.
  }
  cout << nAccH << mAccH;

return 0;
}
