// test911.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Particle Physics Phenomenology 2015 task 4.
// Study jet clustering in kT, C/A and anti-kT.
// Pythia is only used for utilities.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Distance parameter of clustering.
  double Rsep = 0.6;

  // Loop over three clustering algorithms.
  for (int iClus = 0; iClus < 3; ++iClus) {
    cout << "\n\n clustering method = " << iClus << endl;
    if (iClus == 0) cout << fixed << setprecision(1);
    if (iClus == 1) cout << fixed << setprecision(3);
    if (iClus == 2) cout << fixed << setprecision(5);

    // Information arrays and values.
    // status = 1: exists; = 2: joined; = 3: final jet.
    int status[10];
    double y[10], pT[10], pZ[10], E[10];
    int i1Min, i2Min;
    double dij, dijMin;

    // Initial configuration.
    status[1] = 1; y[1] = -0.5; pT[1] = 10.;
    status[2] = 1; y[2] =  0.0; pT[2] = 80.;
    status[3] = 1; y[3] =  0.4; pT[3] = 60.;
    //status[4] = 1; y[4] =  0.7; pT[4] = 20.;
    int nPart = 3;
    for (int i = 1; i <= nPart; ++i) {
      pZ[i] = pT[i] * sinh(y[i]);
      E[i]  = pT[i] * cosh(y[i]);
    }

    // Loop over clustering steps.
    for ( ; ; ) {
      i1Min  = 0;
      i2Min  = 0;
      dijMin = 1e10;

      // Loop over particles. Distance to beams.
      for ( int i1 = 1; i1 <= nPart; ++i1)
      if (status[i1] == 1) {
        if      (iClus == 0) dij = pT[i1];
        else if (iClus == 1) dij = 1;
        else                 dij = 1./pT[i1];
        if (dij < dijMin) {
          i1Min  = i1;
          i2Min  = 0;
          dijMin = dij;
        }
        cout << " i1 = " << i1 << " i2 = 0 dij = " << dij << endl;

        // Distance between particles.
        for ( int i2 = i1 + 1; i2 <= nPart; ++i2)
        if (status[i2] == 1) {
          if      (iClus == 0) dij = min( pT[i1], pT[i2]);
          else if (iClus == 1) dij = 1.;
          else                 dij = min( 1./pT[i1], 1./pT[i2]);
          dij *= abs(y[i1] - y[i2]) / Rsep;
          if (dij < dijMin) {
            i1Min  = i1;
            i2Min  = i2;
            dijMin = dij;
          }
          cout << " i1 = " << i1 << " i2 = " << i2 << " dij = " << dij
               << endl;

        // End loop over particles and pairs.
        }
      }

      // Done if nothing found.
      if (i1Min == 0) break;

      // Remove jet from list.
      if (i2Min == 0) {
        status[i1Min] = 3;
        cout << " >>> Jet " << i1Min << " removed from list" << endl;

      // Combine two jets into new one.
      } else {
        ++nPart;
        status[i1Min] = 2;
        status[i2Min] = 2;
        status[nPart] = 1;
        pT[nPart] = pT[i1Min] + pT[i2Min];
        pZ[nPart] = pZ[i1Min] + pZ[i2Min];
        E[nPart]  =  E[i1Min] +  E[i2Min];
        y[nPart]  = 0.5 * log( (E[nPart] + pZ[nPart]) /
                    (E[nPart] - pZ[nPart]) );
        cout << " >>> Jets " << i1Min << " and " << i2Min
             << " combined to new jet " <<  nPart << endl;
      }

    // End loop over clustering steps.
    }

    // List of clustering.
    for (int i = 1; i <= nPart; ++i) cout << " i = " <<  i << " status "
      << status[i] << " y = "  << setprecision(3) << y[i]
      << setprecision(1) << " pT = " << pT[i]  << " pZ = "  << pZ[i]
      <<  " E = "  << E[i] << endl;


  // End loop over three clustering algorithms.
  }

  // Done.
  return 0;
}
