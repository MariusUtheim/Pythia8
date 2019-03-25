// test185.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Generate a second hard interaction by MPI machinery, to check
// if enhancement factor catches correct behaviour. Cf. test181.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  //Event& event = pythia.event;

  // Set up hard process.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 30.");

  // Set up MPIs.
  pythia.readString("MultipartonInteractions:bProfile = 3");
  pythia.readString("MultipartonInteractions:expPow = 1.");
  pythia.readString("MultipartonInteractions:processLevel = 1");

  // Switch off other things.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:Remnants = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Check:event = off");
  pythia.readString("Next:numberCount = 100000");

  // Initialize for LHC at 13 TeV.
  pythia.readString("Beams:eCM = 13000.");
  pythia.init();

  // Histograms.
  Hist pTfirst("pT first collision",    100, 0., 400.);
  Hist pTsecond("pT second collision",  100, 0., 200.);
  Hist nMult("number of multiparton interactions", 100, -0.5, 99.5);
  Hist bMore("b enhancement factor",    100, 0., 20.);

  // Statistics.
  int    nev = 1000000;
  int    n2i = 0;
  int    n50 = 0;
  double eal = 0.;
  double e2i = 0.;
  double e50 = 0.;

  // Generate events.
  for (int iev = 0; iev < nev; ++iev) {
    pythia.next();

    // Histogram multiparton interactions.
    int    nMPI = pythia.info.nMPI();
    double eMPI = pythia.info.enhanceMPI();
    nMult.fill( nMPI );
    bMore.fill( eMPI);
    eal += eMPI;

    // Histogram first and second pT. Statistics pT > 50.
    double pT1 = pythia.info.pTMPI(0);
    pTfirst.fill( pT1 );
    if (nMPI > 1) {
      double pT2 = pythia.info.pTMPI(1);
      pTsecond.fill( pT2 );
      ++n2i;
      e2i += eMPI;
      if (pT2 > 30.) {
        ++n50;
        e50 += eMPI;
      }
    }

  }

  // Statistics.
  pythia.stat();
  cout << " number of events  = " << n50 << endl;
  double sigma1  = pythia.info.sigmaGen();
  double sigmaND = 56.79;
  double sigma2  = pow2(sigma1) / (2. * sigmaND);
  double nNew    = nev * sigma2 / sigma1;
  double enhal   = eal / max( 1, nev);
  double enh2i   = e2i / max( 1, n2i);
  double enh50   = e50 / max( 1, n50);
  cout << fixed << setprecision(3) << " enhal = " << enhal << " enh2i = "
       << enh2i << " enh50 = " << enh50 << endl;
  cout << " analytical enhancement = " << pythia.info.enhanceMPIavg() << endl;
  cout << fixed << setprecision(2) << " expected without enhancement = "
       << nNew << endl;
  cout << fixed << setprecision(2) << " expected with al enhancement = "
       << nNew * enhal << endl;
  cout << fixed << setprecision(2) << " expected with 2i enhancement = "
       << nNew * enh2i << endl;
  cout << fixed << setprecision(2) << " expected with 50 enhancement = "
       << nNew * enh50 << endl;

  // Print histograms.
  cout << pTfirst << pTsecond << nMult << bMore;

  // Done.
  return 0;
}
