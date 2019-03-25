// test238.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how to extract the PDF of a second hard interaction,
// given the nature of the first interaction. Check on momentum sum rule.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Integration over the second PDF for given first PDF.
  // Number of points, x ranges and local variables.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double x2, sumNow;

  // Generator needs to be initialized, for some dummy process,
  // to have everything initialized. Also place to pick PDF set.
  Pythia pythia;
  pythia.readString("HardQCD:hardbbbar = on");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PDF:pSet = 8");
  pythia.init();

  // Generate one event to set up structure.
  //pythia.next();

  // Extract shorthand for the first beam particle.
  BeamParticle& beam = pythia.beamA;

  // Set Q2 scales for the two interactions.
  double Q12 = 100.;
  double Q22 = 20.;

  // Loop over values for the first interaction.
  for (int iid1 = 1; iid1 < 6; ++iid1)
  //for (int iid1 = -5; iid1 < 6; ++iid1)
  for (int ix1 = 0; ix1 < 7; ++ ix1) {
    int id1    = (iid1 == 0) ? 21 : iid1;
    double x1  = (ix1 < 4) ? pow( 10., ix1 - 4) : 0.1 * (ix1 - 2);

    // Set up beams for first interaction.
    // For u and d quarks the valence or sea assignment for the first
    // interaction is random based on relative PDF at given x and Q2.
    // Therefore redo pickValSeaComp inside loop for right average.
    beam.clear();
    beam.append( 3, id1, x1);
    beam.xfISR(  0, id1, x1, Q12);
    beam.pickValSeaComp();

    // Initialize momentum sum to begin integration.
    double sum = 0.;

    // Integration at large x in linear steps.
    for (int iLin = 0; iLin < nLin; ++iLin) {
      if (id1 == 1 || id1 == 2) beam.pickValSeaComp();
      x2     = xLin + (iLin + 0.5) * dxLin;
      sumNow = beam.xfMPI( 21, x2, Q22) + beam.xfMPI( 22, x2, Q22);
      for (int i = 1; i < 6; ++i)
        sumNow += beam.xfMPI( i, x2, Q22) + beam.xfMPI( -i, x2, Q22);
      sum   += dxLin * sumNow;
    }

    // Integration at small x in logarithmic steps.
    for (int iLog = 0; iLog < nLog; ++iLog) {
      if (id1 == 1 || id1 == 2) beam.pickValSeaComp();
      x2     = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
      sumNow = beam.xfMPI( 21, x2, Q22) + beam.xfMPI( 22, x2, Q22);
      for (int i = 1; i < 6; ++i)
        sumNow += beam.xfMPI( i, x2, Q22) + beam.xfMPI( -i, x2, Q22);
      sum   += dxLog * x2 * sumNow;
    }

    // Result.
    cout << " For id1 = " << setw(2) << id1 << " and x1 = " << fixed
         << setprecision(5) << x1 << " the momentum sum is "
         << setprecision(5) << x1 + sum << endl;
  }

  // Done.
  return 0;
}
