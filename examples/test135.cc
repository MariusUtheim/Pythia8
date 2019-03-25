// main29.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Omnibus version for colour reconnection (CR) effect studies.

// Links to two UserHooks that, along with the internal models,
// implement all the models used for the top mass study in
// S. Argyropoulos and T. Sjostrand,
// arXiv:1407.6653 [hep-ph] (LU TP 14-23, DESY 14-134, MCnet-14-15)

// Warning: some small modifications have been made when collecting
// the models, but nothing intended to change the behaviour.
// Note: the move model is also available with ColourReconnection:mode = 2,
// while the ColourReconnection:mode = 1 model has not been used here.
// Note: the new models tend to be slower than the default CR scenario,
// since they have to probe many more reconnection possibilities.

// Important: the top mass shift analysis encoded here is very primitive,
// does not perform well at all, and should not be taken seriously.
// The important part is that you see how the different scenarios
// should be set up to operate as intended.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate.
  // Warning: much statistics is needed for significant results,
  // so this is just an appetizer. Anyway, the reconstruction is
  // pretty lousy, so not useful for real studies.
  int nEvent = 1000;

    // Generator at 8 TeV LHC.
    Pythia pythia;
    //Event& event = pythia.event;
    pythia.readString("Beams:eCM = 8000.");

    // g g -> t tbar.
    pythia.readString("Top:gg2ttbar = on");

    // Simplify generation. For tryout only.
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:FSR = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("BeamRemnants:primordialKT = off");
    pythia.readString("HadronLevel:all = off");

    // Top and W masses. Semileptonic top decay chosen by W decay.
    // (One of two charge states, so properly ought to symmetrize.)
    // Trick: only allow decay to stable tau, standing in for e and mu
    // as well, but the tau is easy to remove before jet finding.
    pythia.readString("6:m0 = 173.3");
    pythia.readString("24:m0 = 80.385");
    pythia.readString("24:onMode = off");

    // This setup still produces only muonic W- decays.
    pythia.readString("24:onIfAny = 11");
    pythia.readString("24:onNegIfAny = 13");

    // This setup still produces only muonic W- decays.
    //pythia.readString("24:onNegIfAny = 11 12 13 14");
    //pythia.readString("24:onPosIfAny = 11 12");

    // This setup works, but is not quite the desired result.
    //pythia.readString("24:onNegIfAny = 13 14");
    //pythia.readString("24:onPosIfAny = 11 12");

    // This setup crashes.
    //pythia.readString("24:onNegIfAny = 11 12");
    //pythia.readString("24:onPosIfAny = 11 12");

    // Reduce printout.
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberCount = 100000");

    // Initialize.
    pythia.init();

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      pythia.event.list();

      // End of event loop. Statistics. Histograms.
    }

    pythia.stat();

  // Done.
  return 0;
}
