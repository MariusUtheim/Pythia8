// test235.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("test235.cmnd");

  // Extract settings to be used in the main program.
  int nEvent    = pythia.mode("Main:numberOfEvents");
  int nAbort    = pythia.mode("Main:timesAllowErrors");
  int idChi     = 1000022;

  // Initialization.
  pythia.init();

  // Histograms.
  Hist chiPro( "chi production distance ", 90, -6., 3.);
  Hist chiDec( "chi decay distance ", 90, -6., 3.);
  Hist chiDau( "chi daughter production distance ", 90, -6., 3.);
  Hist vtxMis( "relative mismatch chi vs. primary hadrons", 100, -12., -2.);

  // Begin event loop.
  int iAbort    = 0;
  int nErrAtlas = 0;
  int nMatch    = 0;
  int nMismatch = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    //cout<<" EVENT "<<iEvent<<endl;

    // Generate event.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
        cout << " Aborted since reached end of Les Houches Event File\n";
        break;
      }

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Histogram production and decay vertices of neutralino.
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].id() != idChi) continue;
      int dau1 = event[i].daughter1();
      if (dau1 == 0) { cout << " Error: no daughter " << endl; continue;}
      if (event[dau1].id() == idChi) continue;
      double cPro = log10( max(1e-7, event[i].vProd().pAbs() ) );
      double cDec = log10( max(1e-7, event[i].vDec().pAbs() ) );
      double cDau = log10( max(1e-7, event[dau1].vProd().pAbs() ) );
      chiPro.fill( cPro);
      chiDec.fill( cDec);
      chiDau.fill( cDau);

      // Error printout if very short-lived neutralino.
      Vec4 vDec = event[i].vDec();
      if (vDec.pAbs()  < 1e-6) cout << " in event # " << iEvent
        << scientific << " distance = " << vDec.pAbs()
        << " lifetime = " << event[i].tau() << endl;

      // Error printout as in ATLAS Python code.
      if (abs(vDec.px()) < 0.01 && abs(vDec.py()) < 0.01
        && abs(vDec.pz()) < 0.01) {
        ++nErrAtlas;
       //cout << " in event # " << iEvent << scientific << " x = " << vDec.px()
       //     << " y = " << vDec.py() << " z = " << vDec.pz() << endl;
      }
    }

    // Check that primary hadrons from neutralinos have consistent vertices.
    //bool hasMismatch = false;
    double maxDiff = 0.;
    int iMaxDiff = 0;
    for (int i = 3; i < event.size(); ++i) if (event[i].statusAbs()/10 == 8) {
      int mot1 = i;
      do { mot1 = event[mot1].mother1(); }
      while (mot1 > 0 && mot1 < event.size() && event[mot1].id() != idChi
        && event[mot1].statusAbs()/10 != 9);
      if (event[mot1].id() == idChi) {
        Vec4 vHad = event[i].vProd();
        Vec4 vChi = event[mot1].vDec();
        Vec4 vDif = vHad - vChi;
        double diff = (abs(vDif.px()) + abs(vDif.py()) + abs(vDif.pz())
          + abs(vDif.e()) ) / vChi.e();
        if (diff < 1e-5) ++nMatch;
        else {
          ++nMismatch;
          //if (!hasMismatch) cout << " particle = " << mot1 << " has xDec = "
          // << scientific << vChi.px() << " yDec = " << vChi.py()
          // << " zDec = " << vChi.pz()   << " tDec = " << vChi.e() << endl;
          // hasMismatch = true;
          //cout << " problem in event = " << iEvent << " particle = "
          //  << i << " error = " << diff << endl;
        }
        if (diff > maxDiff) {maxDiff = diff; iMaxDiff = i;}
      }
    }
    //if (hasMismatch) event.list(true);
    vtxMis.fill( max( -11., log10(maxDiff) ) );
    if (maxDiff > 1e-2) {
      cout << " Large maxDiff = " << maxDiff << "for particle "
           << iMaxDiff << endl;
      pythia.process.list(true);
      pythia.event.list(true);
   }

  // End of event loop. Statistics.
  }
  pythia.stat();
  cout << chiPro << chiDec << chiDau << vtxMis;
  cout << "\n Number of non-displaced neutralino decays (ATLAS definition) = "
       << nErrAtlas << endl;
  cout << " Number of primary hadrons with matching production vertices = "
       << nMatch << "\n and with mismatching ones = " << nMismatch << endl;

  // Done.
  return 0;
}
