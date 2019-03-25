// test187.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check method(s) to find transverse vectors.
// Outcome: my method slightly worse accuracy and factor ~5 slower.
// Create two vectors that are perpendicular to both input vectors.

//==========================================================================

// in Basics.h:

  // Create two vectors that are perpendicular to both input vectors.
  friend pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2,
    bool noBoost);

// Create two vectors that are perpendicular to both input vectors.
pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2,
  bool noBoost = false);

// in Basics.cc:

pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2,
bool noBoost) {

  if (noBoost) {
    // One perpendicular vector from three-dimensional cross-product.
    Vec4 nPerp( cross3(v1,v2) );
    double TINY = std::numeric_limits<double>::epsilon();
    if ( abs(nPerp.pAbs()) < TINY) {
      Vec4 aux;
      if (v1.px() != 0.)      aux.p(v1.yy,v1.px(),v1.pz(),v1.e());
      else if (v1.py() != 0.) aux.p(v1.px(),v1.pz(),v1.py(),v1.e());
      else if (v1.pz() != 0.) aux.p(v1.pz(),v1.py(),v1.px(),v1.e());
      nPerp.p( cross3(v1,aux) );
    }
    nPerp /= abs(nPerp.pAbs());

    // Second perpendicular vector from four-dimensional cross-product.
    Vec4 lPerp( cross4(v1,v2,nPerp) );
    lPerp /= sqrt(abs(lPerp.m2Calc()));
    return make_pair(nPerp,lPerp);

  // Alternative method that uses boost from rest frame.
  } else {
    Vec4 eX( 1., 0., 0., 0.);
    Vec4 eY( 0., 1., 0., 0.);
    RotBstMatrix tov1v2;
    tov1v2.fromCMframe( v1, v2);
    eX.rotbst(tov1v2);
    eY.rotbst(tov1v2);
    return make_pair( eX, eY );
  }
}

//==========================================================================

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Extract random number generator. Histograms.
  Pythia pythia;
  Rndm rndm = pythia.rndm;
  Hist errH1( "Error method 1", 100, 0., 1e-14);
  Hist errH2( "Error method 2", 100, 0., 1e-14);

  // Fix configuration.
  /*
  Vec4 v1, v2;
  do {
    v1 = Vec4( rndm.flat(), rndm.flat(), rndm.flat(), 3. * rndm.flat() );
    v2 = Vec4( rndm.flat(), rndm.flat(), rndm.flat(), 3. * rndm.flat() );
  } while ((v1 + v2).m2Calc() < 0.);
  */

  // Loop over configurations.
  for (int iEvent = 0; iEvent < 100000000; ++iEvent) {

    // Two random vectors.
    Vec4 v1( rndm.flat(), rndm.flat(), rndm.flat(), 3. * rndm.flat() );
    Vec4 v2( rndm.flat(), rndm.flat(), rndm.flat(), 3. * rndm.flat() );
    if ( (v1 + v2).m2Calc() < 0.) continue;

    // Stefan's method.
    pair<Vec4,Vec4> method1 = getTwoPerpendicular( v1, v2);
    Vec4 eX1 = method1.first;
    Vec4 eY1 = method1.second;

    // My method.
    pair<Vec4,Vec4> method2 = getTwoPerpendicular( v1, v2, false);
    Vec4 eX2 = method2.first;
    Vec4 eY2 = method2.second;

    // Compare outcome: four-product checks.
    double err1  = max( max( max(abs(v1*eX1), abs(v1*eY1)),
                             max(abs(v2*eX1), abs(v2*eY1)) ),
                        max( max(abs(eX1*eX1+1.), abs(eY1*eY1+1.)),
                                 abs(eX1*eY1) ) ) ;
    double err2  = max( max( max(abs(v1*eX2), abs(v1*eY2)),
                             max(abs(v2*eX2), abs(v2*eY2)) ),
                        max( max(abs(eX2*eX2+1.), abs(eY2*eY2+1.)),
                                 abs(eX2*eY2) ) ) ;
    errH1.fill( err1 );
    errH2.fill( err2 );

    // Compare outcome: vector and four-product printout.
    if (iEvent < 10) {
      cout << "\n  v1 = " << v1 << "  v2 = " << v2
           << " eX1 = " << eX1 << " eY1 = " << eY1
           << " eX2 = " << eX2 << " eY2 = " << eY2;
      cout << scientific << setprecision(5) << " err1 = " << err1
           << " err2 = " << err2 << endl;
    }
  }

  // Print histograms and done.
  cout << errH1 << errH2;
  return 0;
}
