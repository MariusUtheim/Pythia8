// test909.cc
// Violations of total cross sections.
// Exercise 9.2.

#include "../include/Pythia8/PythiaStdlib.h"
using namespace Pythia8;

int main() {

  // Conversion factor and header.
  double conv = 0.0511;
  cout << "     e          s       sigTot      Bel1      sigEl1      "
       << "Bel2      sigEl2" << endl;

  // Loop over different energies.
  for (int i = 1; i < 50; ++i) {
    double e = pow( 10., 0.5*i);
    double s = e * e;

    // Cross sections and slopes.
    double sigmatot = 21.7 * pow(s, 0.08);
    double Bel1     = 9.2 + 0.5 * log( 0.25 * s);
    double Bel2     = 12.0 - 0.22 * log(s) + 0.037 * log(s) * log(s);
    double sigmael1 = conv * sigmatot * sigmatot / Bel1;
    double sigmael2 = conv * sigmatot * sigmatot / Bel2;

    // Printout.
    cout << scientific << setprecision(3) << " " << e << "  " << s
         << "  " << sigmatot << "  " << Bel1 << "  " << sigmael1
         << "  " << Bel2 << "  " << sigmael2 << endl;

  // Done.
  }
  return 0;
}
