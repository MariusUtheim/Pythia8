// mymain02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program for the pT_Z spectrum at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Z with 0 or 1 jet. Minimal pT for latter. CM energy, number of events.
  int nJet     = 1;
  double pTmin = 20.;
  double eCM   = 13000.;
  int nEvent   = 10000;

  // Generator.
  Pythia pythia;

  // Set CM energy of collisions.
  pythia.settings.parm("Beams:eCM", eCM);

  // Select process and set pTmin.
  if (nJet == 0) {
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  } else {
    pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
    pythia.readString("WeakBosonAndParton:qg2gmZq = on");
    pythia.settings.parm("PhaseSpace:pTHatMin", pTmin);
  }

  // Restrict range of Z0 mass. Force decay to electron pair.
  pythia.readString("23:mMin = 76.");
  pythia.readString("23:mMax = 106.");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 11");

  // Switch off everything except the hard process and beam remnants.
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("BeamRemnants:primordialKT = off");
  pythia.readString("HadronLevel:all = off");

  // Initialize. Sum of Sudakov weights.
  pythia.init();
  double sudSum = 0.;

  // Histogram pT values.
  Hist pTZnoSud("dN/dpTZ without Sudakov", 100, 0., 100.);
  Hist pTZwithSud("dN/dpTZ with Sudakov", 100, 0., 100.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Loop over particles in event. Find last Z0 copy and its pT.
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 23) iZ = i;
    double pTZ = pythia.event[iZ].pT();

    // Use the Info methods to read off main process characteristics.
    double pTnow = pythia.info.pTHat();
    double x1now = pythia.info.x1();
    double x2now = pythia.info.x2();
    int id1now   = pythia.info.id1();
    int id2now   = pythia.info.id2();

    // Check that consistent results for pTZ.
    if (abs(pTnow - pTZ) > 0.1) cout << " Mismatch in pT of Z, pTZ = "
      << fixed << setprecision(3) << pTZ << " vs. pTnow = " << pTnow << endl;

    // Calculate Sudakov factors. Purely fictitious for now.
    double Q2   = pTnow * pTnow;
    double sud1 = exp(-0.1 * log(0.01 * eCM / pTnow)) * (1. - x1now);
    if (id1now == 21) sud1 *= sud1;
    double sud2 = exp(-0.1 * log(0.01 * eCM / pTnow)) * (1. - x2now);
    if (id2now == 21) sud2 *= sud2;
    double sud = sud1 * sud2;
    sudSum += sud;

    // Histogram pT spectrum without or with Sudakov.
    pTZnoSud.fill( pTnow );
    pTZwithSud.fill( pTnow, sud );

  // End of event loop. Statistics, especially cross section (in pb).
  }
  pythia.stat();
  double sigmaTot = pythia.info.sigmaGen() * 1e9;
  double sigmaSud = sigmaTot * sudSum / nEvent;
  cout << endl << fixed << setprecision(3)
       << " Total cross section = " << sigmaTot << " pb" << endl
       << " Ditto, with Sudakov = " << sigmaSud << " pb" << endl;

  // Simple and somewhat better histogram output.
  cout << pTZnoSud << pTZwithSud;
  HistPlot hpl("mymain02");
  hpl.frame( "pTZ", "Z0 pT spectrum", "pT (GeV)", "dN/dpTZ");
  hpl.add( pTZnoSud, ",b");
  hpl.add( pTZwithSud, ",r");
  hpl.plot();

  // Done.
  return 0;
}
