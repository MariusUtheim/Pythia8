// test916.cc
// Integration of MPI rate to get cross section.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//-----------------------------------------------------------------

int main() {

  // Set parameter values.
  double bMax    = 10.;
  int    nStep   = 100000;
  double b0      = 1.;
  double nMax[5] = {1., 2., 5., 10., 20.};
  double sum[5]  = {0., 0., 0., 0., 0.};
  double pNow;

  // Loop over integration range and  evaluate inner exponential.
  for (int i = 0; i < nStep; ++i) {
    double bNow = bMax * (i + 0.5) / nStep;
    double bExp = exp(- pow2(bNow / b0));

    // Loop over nMax values and evaluate contribution to integral.
    for (int iMax = 0; iMax < 5; ++iMax) {
      pNow       = 1 - exp(- nMax[iMax] * bExp);
      sum[iMax] += bNow * pNow;
    }
  }

  // Normalize integrals and print values.
  double norm = 2. * M_PI * bMax / nStep;
  for (int iMax = 0; iMax < 5; ++iMax)
    cout << " nMax = " << fixed << setprecision(2) << nMax[iMax]
         << " sum = " << setprecision(3) << norm * sum[iMax] << endl;

  // Done.
  return 0;
}
