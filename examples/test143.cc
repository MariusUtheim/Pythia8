// test143:
// Check equivalence of matrix elements in WeakShowerMEs.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME uG -> uG, with an overall factor of
// g_s^4 / 9 was removed.

double getMEqg2qg(double sH, double tH, double uH, bool useNew) {

  double sH2 = sH* sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  return (sH2 + uH2) * (9. / tH2 - 4. / (sH * uH));

  } else {
  double sH3 = sH2 * sH;
  return (18.*sH3*uH - 4*tH2*uH2 + 9*sH*uH2*(tH + 2*uH) +
    sH2*(-4.*tH2 + 9*tH*uH + 18*uH2))/(sH*tH2*uH);
  }
}

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME ud -> ud, with an overall factor of
// g_s^4 / 9 was removed.

double getMEqq2qq(double sH, double tH, double uH, bool sameID, bool useNew) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  if (sameID) return 2. * ( (sH2 + uH2) / tH2 + (sH2 + tH2) / uH2
    - 2. * sH2 / (3. * tH * uH) );
  else return 4. * (sH2 + uH2) / tH2;

  } else {

  if (sameID) return 4./2.*((sH2+uH2)/tH2 + (sH2+tH2)/uH2 -
                            2.*sH2/(3.*tH*uH));
  else return 4.*(sH2+uH2)/tH2;
  }
}

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME gg -> gg, with an overall factor of
// g_s^4 / 9 was removed.

double getMEgg2gg(double sH, double tH, double uH, bool useNew) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  return (81. / 8.) * ( (tH2 + uH2) / sH2 + (sH2 + uH2) / tH2
    + (sH2 + tH2) / uH2 + 3. );

  } else {
  double sigTS  = (9./4.) * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH
                  + sH2 / tH2);
  double sigUS  = (9./4.) * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH
                  + sH2 / uH2);
  double sigTU  = (9./4.) * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH
                  + uH2 / tH2);
  // 0.5 from identical gluons and 9 to get common factor.
  return  0.5 * 9 * (sigTS + sigUS + sigTU);
  }
}

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME gg -> qqbar, with an overall factor of
// g_s^4 / 9 was removed.

double getMEgg2qqbar(double sH, double tH, double uH, bool useNew) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  return (tH2 + uH2) * ( 3. / (2. * tH * uH) - 27. / (8. * sH2) );

  } else {
  double sigTS = (1./6.) * uH / tH - (3./8.) * uH2 / sH2;
  double sigUS = (1./6.) * tH / uH - (3./8.) * tH2 / sH2;
  return 9. * (sigTS + sigUS);
  }
}

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME qqbar -> qqbar, with an overall factor of
// g_s^4 / 9 was removed.

double getMEqqbar2gg(double sH, double tH, double uH, bool useNew) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  return (tH2 + uH2) * ( 16. / (3. * tH * uH) - 12. / sH2 );

  } else {
  double sigTS  = (32./3.) * uH / tH - (24.) * uH2 / sH2;
  double sigUS  = (32./3.) * tH / uH - (24.) * tH2 / sH2;
  // 0.5 from identical gluons
  return 0.5 * (sigTS + sigUS);
  }
}

//--------------------------------------------------------------------------

// Calculate the 2 to 2 ME gg -> qqbar, with an overall factor of
// g_s^4 / 9 was removed.

double getMEqqbar2qqbar(double sH, double tH, double uH, bool sameID,
  bool useNew) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  if (useNew) {
  if (sameID) return 4. * (tH2 + uH2) / sH2 - (8./3.) * uH2 / (sH * tH)
    + 4. * (sH2 + uH2) / tH2;
  else return 4. * (tH2 + uH2) / sH2;

  } else {
  double sigT   = (4.) * (sH2 + uH2) / tH2;
  double sigST  = - (8./3.) * uH2 / (sH * tH);
  if (sameID) return sigT + sigST;
  else return 4. * (tH2 + uH2) / sH2;
  }

}

//--------------------------------------------------------------------------

int main() {

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  // Ten different phase space points.
  for (int iPS = 0; iPS < 10; ++iPS) {
    double sH = 10. * pythia.rndm.flat();
    double tH = -sH * pythia.rndm.flat();
    double uH = - (sH + tH);
    cout << scientific << setprecision(5) << "\n sH = " << sH
         << " tH = " << tH << " uH = " << uH << endl;
    cout << " A  " << getMEqg2qg( sH, tH, uH, true)
         << " vs " << getMEqg2qg( sH, tH, uH, false) << endl;
    cout << " B1 " << getMEqq2qq( sH, tH, uH, true, true)
         << " vs " << getMEqq2qq( sH, tH, uH, true, false) << endl;
    cout << " B2 " << getMEqq2qq( sH, tH, uH, false, true)
         << " vs " << getMEqq2qq( sH, tH, uH, false, false) << endl;
    cout << " C  " << getMEgg2gg( sH, tH, uH, true)
         << " vs " << getMEgg2gg( sH, tH, uH, false) << endl;
    cout << " D  " << getMEgg2qqbar( sH, tH, uH, true)
         << " vs " << getMEgg2qqbar( sH, tH, uH, false) << endl;
    cout << " E  " << getMEqqbar2gg( sH, tH, uH, true)
         << " vs " << getMEqqbar2gg( sH, tH, uH, false) << endl;
    cout << " F1 " << getMEqqbar2qqbar( sH, tH, uH, true, true)
         << " vs " << getMEqqbar2qqbar( sH, tH, uH, true, false)
         << " note missing s-channel contribution" << endl;
    cout << " F2 " << getMEqqbar2qqbar( sH, tH, uH, false, true)
         << " vs " << getMEqqbar2qqbar( sH, tH, uH, false, false) << endl;
  }

  return 0;
}
