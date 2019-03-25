// test906.cc
// Evolve a quark by ISR.
// Exercise 5.1.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Set up Pythia for use as random number generator.
  Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 32133");
  pythia.init();

  // Loop over three zMax values.
  for (int izMax = 0; izMax < 5; ++izMax) {

  // Values for simulation.
  int nEvent = 10000000;
  double Qmax  = 100.;
  double Qmin = 1.;
  double zMax = 0.9;
  if (izMax == 1) zMax = 0.99;
  if (izMax == 2) zMax = 0.999;
  if (izMax == 3) zMax = 0.9999;
  if (izMax == 4) zMax = 0.99999;
  double alphas = 0.15;
  double xBeg = 0.5;

  // Upper estimate of branching kernel gives z integral.
  double intZ = (8./3.) * (-log(1. - zMax));
  double evol = M_PI / (alphas * intZ);

  // Statistics;
  double xAvg = 0.;
  Hist xDist( "x value after evolution", 100, 0., 1.0);

  // Loop over "events". Initial values.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    double xNow = xBeg;
    double Qnow = Qmin;

    // Iterate upwards in Q.
    for ( ;  ; ) {
      Qnow /= pow(pythia.rndm.flat(), evol);
      if (Qnow > Qmax) break;

      // Pick a z and reject down to correct kernel.
      double z = 1. - pow(1. - zMax, pythia.rndm.flat());
      if (1 + pow2(z) < 2. * pythia.rndm.flat()) continue;

      // Aceptable branching.
      xNow *= z;

    // End of branching and event loops.
    }
    xAvg += xNow;
    xDist.fill( xNow);
  }

  // Statistics.
  cout << xDist;
  cout << " On average the final x value is " << fixed << setprecision(5)
       << xAvg/double(nEvent) << " for zMax = " << zMax << endl;

  // End loop over three zMax values.
  }

  return 0;
}
