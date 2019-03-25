// test226.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Rate of forward protons, for Valery Khoze.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Phase space cuts. Event types.
  int    nEvent  = 1000000;
  double xiMin   = 0.03;
  double xiMax   = 0.15;
  double pTmin   = 0.4;
  double etaMax  = 2.5;
  string types[5] = { "nondiff", "elastic", "single diff", "double diff",
                      "central diff"};

  // Loop over cases. Generator. Reset statistics.
  for (int ica = 0; ica < 2; ++ica) {
    Pythia pythia;
    Event& event = pythia.event;
    int countAll[5] = { 0, 0, 0, 0, 0};
    int countAcc[5] = { 0, 0, 0, 0, 0};

    // Common settings.
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("SoftQCD:all = on");
    pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.readString("Next:numberShowEvent = 0;");
    pythia.readString("Next:numberCount = 10000");

    // Case-by-case settings. 0: SaS; 1: ABMST.
    if (ica == 0) pythia.readString("SigmaTotal:mode = 1");
    if (ica == 1) pythia.readString("SigmaTotal:mode = 3");
    if (ica == 0) pythia.readString("SigmaDiffractive:mode = 1");
    if (ica == 1) pythia.readString("SigmaDiffractive:mode = 3");

    // Initialize.
    pythia.init();

    // Begin event loop.
    int nAbort = 5;
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Quit if too many failures.
      if (!pythia.next()) {
        if (++iAbort < nAbort) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      // Event classification.
      int code1 = pythia.info.code();
      int code2 = 0;
      if      (code1 == 102) code2 = 1;
      else if (code1 == 103 || code1 == 104) code2 = 2;
      else if (code1 == 105) code2 = 3;
      else if (code1 == 106) code2 = 4;

      // Check whether event satisfies criteria.
      bool hascentral = false;
      bool hasp       = false;
      bool has2p      = false;
      for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() ) {
        if (event[i].isCharged() && event[i].pT() > pTmin
          && abs(event[i].eta()) < etaMax) {
          hascentral = true;
          break;
        }
        if (event[i].id() == 2212) {
          double xp = 2. * event[i].e() / event[0].e();
          double xi = 1. - xp;
          if (xi > xiMin && xi < xiMax) {
            if (hasp) has2p = true;
            hasp = true;
          }
        }
      }

      // Statistics on events.
      ++countAll[code2];
      if (!hascentral) {
        if (hasp)  ++countAcc[code2];
        if (has2p) ++countAcc[code2];
      }

    // End of event loop. Print statistics. End case loop.
    }
    pythia.stat();
    double sigmaTot = pythia.info.sigmaGen();
    cout << endl << "          type     total  accepted   sigma (mb)"
      << endl;
    for (int ic = 0; ic < 5; ++ic) cout << setw(14) << types[ic]
      << setw(10) << countAll[ic] << setw(10) << countAcc[ic]
      << setw(12) << fixed << setprecision(5)
      << countAcc[ic] * sigmaTot / nEvent << endl;
  }

  // Done.
  return 0;
}
