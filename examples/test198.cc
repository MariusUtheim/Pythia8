// test198.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Freestanding parametrizations of the SaS cross sections.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// The SigmaSaS class implments the SaS cross sections for pp collisions.

// Coupling constants for DL pp cross sections.
const double EPSILON   = 0.0808;
const double ETA       = -0.4525;
const double XP        = 21.70;
const double YP        = 56.08;
const double BETA0     = 4.658;
const double BHAD      = 2.3;
const double ALPHAPRIME = 0.25;
// Conversion coefficients = 1/(16pi) * (mb <-> GeV^2) * (G_3P)^n,
// with n = 0 elastic, n = 1 single and n = 2 double diffractive.
const double CONVERTEL = 0.0510925;
const double CONVERTSD = 0.0336;
const double CONVERTDD = 0.0084;
const double CRES      = 2.;
const double M2RES     = 4.;
const double SPROTON   = 0.880;

class SigmaSaS {

public:

  // Total cross section.
  double sigmaTot( double s) {
    double sEps = pow( s, EPSILON);
    double sEta = pow( s, ETA);
    return XP * sEps + YP * sEta;
  }

  // Elastic cross section dsigma_el/dt.
  double dsigmaEl( double s, double t) {
    double sigEl0 = CONVERTEL * pow2(sigmaTot(s));
    double bEl    = 4. * BHAD + 4. * pow( s, EPSILON) - 4.2;
    return sigEl0 * exp(bEl * t);
  }

  // Single diffractive cross section.
  double dsigmaSD( double s, double m2, double t) {
    double sigSD0 = CONVERTSD * pow3(BETA0) / m2;
    double fSD    = (1. - m2 / s) * (1. + CRES * M2RES / (M2RES + m2) );
    double bSD    = 2. * BHAD + 2. * ALPHAPRIME * log(s/m2);
    return sigSD0 * fSD * exp(bSD * t);
}

  // Double diffractive cross section.
  double dsigmaDD( double s, double mA2, double mB2, double t) {
    double sigDD0 = CONVERTDD * pow2(BETA0) / (mA2 * mB2);
    double fDD    = (1. - pow2(sqrt(mA2) + sqrt(mB2)) / s)
                  * s * SPROTON / (s * SPROTON + mA2 * mB2)
                  * (1. + CRES * M2RES / (M2RES + mA2) )
                  * (1. + CRES * M2RES / (M2RES + mB2) );
    double bDD    = 2. * ALPHAPRIME
                  * log( exp(4.) + s / (mA2 * mB2 * ALPHAPRIME) );
    return sigDD0 * fDD * exp(bDD * t);
  }

};

//==========================================================================

// The main program.

int main() {

  // Set up class for calculations.
  SigmaSaS sigma;

  // Loop over energies and xi values.
  for (int ie = 2; ie < 6; ++ie)
  for (int ix = 1; ix < 9; ++ix) {

    // Values. Skip if small diffractive mass or no rapidity gap.
    double eCM    = pow( 10., ie);
    double xi     = pow( 0.1, ix);
    double t      = 0.;
    double s      = eCM * eCM;
    double m2     = xi * s;
    if (m2 < 2. || m2 > 0.5 * eCM) continue;

    // Calculate cross sections.
    //double sigTot = sigma.sigmaTot(s);
    double sigEl0 = sigma.dsigmaEl(s, t);
    double sigSD0 = sigma.dsigmaSD(s, m2, t);
    double sigDD0 = sigma.dsigmaDD(s, m2, m2, t);
    //cout << fixed << setprecision(3) << " " << sigTot << " " << sigEl0
    //     << " " << sigSD0 << " " << sigDD0 << endl;

    // Interesting ratio sigma_DD * sigma_el / sigma_SD^2.
    double ratio  = sigDD0 * sigEl0 / pow2(sigSD0);
    double rCorr  = ratio / pow( s / 400., 2. * EPSILON);
    cout << " ratio for eCM = " << scientific << setprecision(2) << eCM
         << " and xi = " << xi << " is = " << fixed << setprecision(3)
         << setw(6) << ratio << " and corrected = " << setw(6) << rCorr
         << endl;
  }

  // Done.
  return 0;

}
