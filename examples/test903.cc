// test903.cc
// Simple shower: a single quark branches consecutively.
// Exercise 3.4

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Set up Pythia for use as random number generator.
  Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 32133");
  pythia.init();

  // Values for simulation.
  int nEvent = 10000;
  double Qmax = 91.;
  double Qmin = 1.;
  double zmax = 0.99;
  double alphas = 0.15;

  // Upper estimate of branching kernel gives z integral.
  double intZ = (8./3.) * (-log(1. - zmax));
  double evol = M_PI / (alphas * intZ);

  // Statistics;
  int nBranch = 0;

  // Loop over "events". Initial values.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    double Qnow = Qmax;

    // Iterate downwards in Q.
    for ( ;  ; ) {
      Qnow *= pow(pythia.rndm.flat(), evol);
      if (Qnow < Qmin) break;

      // Pick a z and reject down to correct kernel.
      double z = 1. - pow(1. - zmax, pythia.rndm.flat());
      if (1 + pow2(z) < 2. * pythia.rndm.flat()) continue;

      // Aceptable branching.
      ++nBranch;

    // End of branching and event loops.
    }
  }

  // Statistics.
  cout << " On average there are " << fixed << setprecision(3)
       << double(nBranch)/double(nEvent) << " branchings per cascade"
       << endl;

  return 0;
}
