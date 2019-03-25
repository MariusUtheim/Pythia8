// test915.cc
// How to set up a random number generator using Pythia.
// Glauber model of nuclear clollisions.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//-----------------------------------------------------------------

// Parameters of the Wood-Saxon distribution (in fm).
int    nNucl  = 208;
double radius = 6.62;
double width  = 0.546;
double rcore  = 0.4;

// Nucleon-nucleon cross section and effectiv max impact parameter.
double sigma  = 6.4;
double dMax   = sqrt(sigma / M_PI);

// Positions of nucleons in the two nuclei and whether wounded.
Vec4 nuclA[208], nuclB[208];
bool woundA[208], woundB[208];

//-----------------------------------------------------------------

// Distribute the nucleons inside a nucleus.

void woodSaxon( int iCase, double bx, double by, Pythia& pythia) {

  // Internal variables.
  double rr, ws, cth, the, phi, xx, yy;

  // Pick one nucleon at a time.
  for (int i = 0; i < nNucl; ++i) {

    // Pick random point inside 2 * radius big sphere.
    do {
      rr = 2. * radius * pow( pythia.rndm.flat(), 1./3.);
      ws = 1. / ( 1. + exp( (rr - radius) / width) );
    } while (ws < pythia.rndm.flat() );
    cth = 2. * pythia.rndm.flat() - 1.;
    the = acos(cth);
    phi = 2. * M_PI * pythia.rndm.flat();

    // Put result in the respective vector.
    xx = bx + rr * sin(the) * cos(phi);
    yy = by + rr * sin(the) * sin(phi);
    if (iCase == 1) nuclA[i] = Vec4( xx, yy, rr * cth, 0.);
    else            nuclB[i] = Vec4( xx, yy, rr * cth, 0.);

  // End of loop over nucleons.
  }

}

//-----------------------------------------------------------------

// Count number of wounded nucleons.

int countWounded() {

  // Internal variables.
  double dNow;

  // Reset nucleon status to not wounded.
  for (int i = 0; i < nNucl; ++i) {
    woundA[i] = false;
    woundB[i] = false;
  }

  // Loop through all pairs of nucleons from side A and side B.
  for (int iA = 0; iA < nNucl; ++iA)
  for (int iB = 0; iB < nNucl; ++iB) {
    dNow = sqrt( pow2(nuclA[iA].px() - nuclB[iB].px())
         +       pow2(nuclA[iA].py() - nuclB[iB].py()) );
    if (dNow < dMax) {
      woundA[iA] = true;
      woundB[iB] = true;
    }
  }

  // Count up number of wounded nucleons and return it.
  int nWound = 0;
  for (int i = 0; i < nNucl; ++i) {
    if (woundA[i]) ++nWound;
    if (woundB[i]) ++nWound;
  }
  return nWound;

}

//-----------------------------------------------------------------

int main() {

  // Set up Pythia for use as random number generator.
  Pythia pythia;
  pythia.init();

  // Histograms.
  Hist nw0( "b =  0.0", 85, 0., 425.1);
  Hist nw1( "b =  2.5", 85, 0., 425.1);
  Hist nw2( "b =  5.0", 85, 0., 425.1);
  Hist nw3( "b =  7.5", 85, 0., 425.1);
  Hist nw4( "b = 10.0", 85, 0., 425.1);
  Hist nw5( "b = 12.5", 85, 0., 425.1);
  Hist nw6( "b = 15.0", 85, 0., 425.1);
  Hist nw7( "0 < b < 20", 85, 0., 425.1);

  // Loop over different impact parameters.
  for (int iImp = 0; iImp < 8; ++iImp) {
    double bImp = 2.5 * iImp;

    // Loop over events.
    int nEvent = 100000;
    double   w = 0.2 / nEvent;
    int nAbove0 = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Select impact parameter. Set up colliding Wood-Saxons.
      if (iImp == 7) bImp = 20. * sqrt(pythia.rndm.flat());
      woodSaxon( 1,  0.5 * bImp, 0., pythia);
      woodSaxon( 2, -0.5 * bImp, 0., pythia);

      // Get number of wounded nucleons. Histogram it.
      int nWound = countWounded();
      if (iImp == 0) nw0.fill( nWound, w);
      if (iImp == 1) nw1.fill( nWound, w);
      if (iImp == 2) nw2.fill( nWound, w);
      if (iImp == 3) nw3.fill( nWound, w);
      if (iImp == 4) nw4.fill( nWound, w);
      if (iImp == 5) nw5.fill( nWound, w);
      if (iImp == 6) nw6.fill( nWound, w);
      if (iImp == 7 && nWound > 0) nw7.fill( nWound, w);
      if (nWound > 0) ++nAbove0;

    // End of event loop and impact parameter loop.
    }
    cout << " iImp = " << iImp << " gives " << nAbove0
         << " nonvanishing" << endl;
  }

  // Show histograms.
  cout << nw0 << nw1 << nw2 << nw3 << nw4 << nw5 << nw6 << nw7;
  HistPlot hpl("test915plot");
  hpl.frame( "out915plot", "Number of wounded nucleons in PbPb collisions",
    "$N_{\\mathrm{wounded}}$", "Probability");
  hpl.add( nw0);
  hpl.add( nw1);
  hpl.add( nw2);
  hpl.add( nw3);
  hpl.add( nw4);
  hpl.add( nw5);
  hpl.add( nw6);
  hpl.plot();
  hpl.frame( "", "Number of wounded nucleons in PbPb collisions",
    "$N_{\\mathrm{wounded}}$", "Probability");
  hpl.add( nw7);
  hpl.plot();
  hpl.frame( "", "Number of wounded nucleons in PbPb collisions",
    "$N_{\\mathrm{wounded}}$", "Probability");
  hpl.add( nw7);
  hpl.plot( true);

  return 0;
}
