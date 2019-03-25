// test901.cc
// How to set up a random number generator using Pythia.
// Exercise 1.3

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Set up Pythia for use as random number generator.
  Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 32133");
  pythia.init();

  // Generate and print a few numbers.
  cout << fixed << setprecision(10);
  for (int i = 0; i < 10; ++i)
    cout << "  " << pythia.rndm.flat() << endl;

  // Study product of three random numbers.
  int nMC = 10000000;
  int nTh = 10000;
  Hist Hmc("distribution of R1 * R2 * R3", 100, 0., 1.);
  Hist Hth("ln^2 x / 2", 100, 0., 1.);
  Hist Hmt("ratio MC/theory", 100, 0., 1.);
  for (int i = 0; i < nMC; ++i) {
    double r3 = pythia.rndm.flat() * pythia.rndm.flat()
              * pythia.rndm.flat();
    Hmc.fill (r3);
  }
  Hmc *= 100. / nMC;
  for (int i = 0; i < nTh; ++i) {
    double x = (i + 0.5) / nTh;
    Hth.fill( x, 0.5 * pow2(log(x)) );
  }
  Hth *= 100. / nTh;
  Hmt = Hmc / Hth;
  cout << Hmc << Hth << Hmt;

  return 0;
}
