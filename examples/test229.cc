// test229.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Proximity of space-time production vertices.
// Output plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Helper class storing info for each final-state hadron.

class HadronRescatterInfo {

public:

HadronRescatterInfo( int iLocIn = 0) : iLoc(iLocIn), nClose(0),
  bClosest(9.9) {}

  int iLoc, nClose;
  double bClosest;

};

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 10000;
  int nAbort  = 5;

  // Maximum impact parameter in fm.
  double sigmaHH = 20.;
  double bMax    = sqrt( 0.1 * sigmaHH / M_PI );

  // Histogram impact parameters and number of possible collisions.
  Hist nHad("number of final hadrons", 100, 0.5, 200.5);
  Hist bClosest("closest impact parameter for hadron", 50, 0., 10.);
  Hist nColl("number of rescatterings per hadron", 10, -0.5, 9.5);

  // Pythia generator. Event record shorthand.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up wanted processes. Vertex selection choices. No decays.
  pythia.readString("Beams:eCM = 2000.");
  pythia.readString("SoftQCD:nondiffractive = on");
  pythia.readString("Fragmentation:setVertices = on");
  //pythia.readString("HadronLevel:Decay = off");
  pythia.readString("111:mayDecay = off");

  // Initialization.
  pythia.init();

  // Array of final hadrons.
  long nTotal = 0;
  vector<HadronRescatterInfo> hadronList;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Loop over all particles and find final hadrons.
    hadronList.clear();
    for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isHadron())
        hadronList.push_back( HadronRescatterInfo(i) );
    int nHadron = hadronList.size();
    nHad.fill( nHadron );
    nTotal += nHadron;

    // Loop over all pairs of final hadrons.
    for (int j1 = 0; j1 < nHadron - 1; ++j1) {
      int i1 = hadronList[j1].iLoc;
      for (int j2 = j1 + 1; j2 < nHadron; ++j2) {
        int i2 = hadronList[j2].iLoc;

        // Extract momentum and vertex info boosted to rest frame of pair.
        Vec4 p1 = event[i1].p();
        Vec4 p2 = event[i2].p();
        Vec4 v1 = event[i1].vProd() * MM2FM;
        Vec4 v2 = event[i2].vProd() * MM2FM;
        RotBstMatrix mRB;
        mRB.toCMframe( p1, p2);
        p1.rotbst(mRB);
        p2.rotbst(mRB);
        v1.rotbst(mRB);
        v2.rotbst(mRB);

        // Impact parameter and motion in z direction.
        double bNow = (v2 - v1).pT();
        double tMax = max( v1.e(), v2.e() );
        double z1   = v1.pz() + (tMax - v1.e()) * p1.pz() / p1.e();
        double z2   = v2.pz() + (tMax - v2.e()) * p2.pz() / p2.e();

        // Store info on impact parameter and number of rescatterings.
        if (z1 < z2) {
          if (bNow < hadronList[j1].bClosest) hadronList[j1].bClosest = bNow;
          if (bNow < hadronList[j2].bClosest) hadronList[j2].bClosest = bNow;
        }
        if (bNow < bMax && z1 < z2) {
          ++hadronList[j1].nClose;
          ++hadronList[j2].nClose;
        }

      // End of pair loops.
      }
    }

    // Histogram impact parameter and nearby hadrons for each hadron.
    for (int j = 0; j < int(hadronList.size()); ++j) {
      bClosest.fill( hadronList[j].bClosest );
      nColl.fill( hadronList[j].nClose );
    }


  // End of event loop. Final statistics.
  }
  pythia.stat();

  // Normalize and print histograms.
  nHad /= nEvent;
  bClosest /= nTotal;
  nColl /= double(nTotal);
  cout << nHad << bClosest << nColl;

  // Write Python code that can generate a PDF file with the distributions.
  HistPlot hpl("test229plot");
  hpl.plotFrame( "out229plot", nHad, "Hadron multiplicity distribution",
    "$n_{\\mathrm{hadrons}}$", "probability");
  hpl.plotFrame( "", bClosest, "Smallest impact parameter",
    "$b_{\\mathrm{closest}}$ (fm)", "probability");
  hpl.plotFrame( "", nColl, "Number of potential rescattering partners",
    "$n_{\\mathrm{partners}}$", "probability");

  // Done.
  return 0;
}
