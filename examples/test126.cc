// test126.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Simple illustration how to implement Dark Matter pair production at the LHC.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// The momentum integral of a companion quark, with its partner at x_s,
// using an approximate gluon density like (1 - x_g)^power / x_g.
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double xCompFrac(int companionPower, double xs) {

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0:
       return xs * ( 5. + xs * (-9. - 2. * xs * (-3. + xs)) + 3. * log(xs) )
         / ( (-1. + xs) * (2. + xs * (-1. + 2. * xs)) );

    case 1:
       return -1. -3. * xs + ( 2. * pow2(-1. + xs) * (1. + xs + xs*xs))
         / ( 2. + xs*xs * (xs - 3.) + 3. * xs * log(xs) );

    case 2:
       return xs * ( (1. - xs) * (19. + xs * (43. + 4. * xs))
         + 6. * log(xs) * (1. + 6. * xs + 4.*xs*xs) ) /
        ( 4. * ( (xs - 1.) * (1. + xs * (4. + xs) )
        - 3. * xs * log(xs) * (1 + xs) ) );

    case 3:
      return 3. * xs * ( (xs - 1.) * (7. + xs * (28. + 13. * xs))
        - 2. * log(xs) * (1. + xs * (9. + 2. * xs * (6. + xs))) )
        / ( 4. + 27. * xs - 31. * pow3(xs)
        + 6. * xs * log(xs) * (3. + 2. * xs * (3.+xs)) );

    default:
      return ( -9. * xs * (xs*xs - 1.) * (5. + xs * (24. + xs)) + 12. * xs
        * log(xs) * (1. + 2. * xs) * (1. + 2. * xs * (5. + 2. * xs)) )
        / ( 8. * (1. + 2. * xs) * ((xs - 1.) * (1. + xs * (10. + xs))
        - 6. * xs * log(xs) * (1. + xs)) );

  }
}

//--------------------------------------------------------------------------

// The x*f pdf of a companion quark at x_c, with its sea partner at x_s,
// using an approximate gluon density like (1 - x_g)^power / x_g.
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double xCompDist(int companionPower, double xc, double xs) {

  // Mother gluon momentum fraction. Check physical limit.
  double xg = xc + xs;
  if (xg > 1.) return 0.;

  // Common factor, including splitting kernel and part of gluon density
  // (and that it is x_c * f that is coded).
  double fac = 3. * xc * xs * (xc*xc + xs*xs) / pow4(xg);

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0:
      return fac / ( 2. - xs * (3. - xs * (3. - 2. * xs)) );

    case 1:
      return fac * (1. - xg) / ( 2. + xs*xs * (-3. + xs) + 3. * xs * log(xs) );

    case 2:
      return fac * pow2(1. - xg) / ( 2. * ((1. - xs) * (1. + xs * (4. + xs))
        + 3. * xs * (1. + xs) * log(xs)) );

    case 3:
      return fac * pow3(1. - xg) * 2. / ( 4. + 27. * xs - 31. * pow3(xs)
        + 6. * xs * log(xs) * (3. + 2. * xs * (3. + xs)) );

    default:
       return fac * pow4(1. - xg) / ( 2. * (1. + 2. * xs) * ((1. - xs)
         * (1. + xs * (10. + xs)) + 6. * xs * log(xs) * (1. + xs)) );

  }
}

//==========================================================================

int main() {

  // Loops.
  for (int cmppw = 0; cmppw < 5; ++cmppw) {
    cout << "\n companion power = " << cmppw << endl
         << scientific << setprecision(8);
    for (int ixs = 1; ixs < 20; ++ixs) {
      double xs = 1. - pow( 0.5, 0.5 * ixs);
      double xCF = xCompFrac( cmppw, xs);
      double xCD = xCompDist( cmppw, 1e-10, xs);
      cout << " xs = " << xs << " xCF = " << xCF << " xCD = " << xCD << endl;
    }
  }

  // Done.
  return 0;
}
