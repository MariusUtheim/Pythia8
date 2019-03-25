// test167.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Radiation off top restricted by top width.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// Global histograms.

Hist egBef( "gluon energies before veto", 100, 0., 10.);
Hist egAft( "gluon energies after veto ", 100, 0., 10.);

//==========================================================================

// Own derived UserHooks class to veto some FSR emissions off top.

class MyUserHooks : public UserHooks {

public:

  // Constructor and destructor do nothing.
  MyUserHooks() {}
  ~MyUserHooks() {}

  // Allow a veto for FSR emissions
  bool canVetoFSREmission() {return true;}

  // Access the event after FSR emission.
  bool doVetoFSREmission(int sizeOld, const Event& event, int , bool ) {

    // Check if emitter is top, else no veto.
    if (event[sizeOld].idAbs() != 6) return false;

    // Invariant mass of t + emitted gluon/photon vs. t itself.
    double mtg    = (event[sizeOld].p() + event[sizeOld + 1].p()).mCalc();
    double deltaM = mtg - event[sizeOld].m();
    egBef.fill( deltaM);

    // Suppression factor for emission or step function.
    double GammaT = particleDataPtr->mWidth(6);
    //double supp = pow2(deltaM) / ( pow2(GammaT) + pow2(deltaM) );
    //if (supp < rndmPtr->flat()) return true;
    if (deltaM < GammaT) return true;
    egAft.fill( deltaM);

    // Else keep event as is.
    return false;
  }

};

//==========================================================================

int main() {

  // Key parameters: # events, cm Energy, t mass, minimum pT.
  int nEvent   = 10000;
  double eCM   = 1000.;
  double mTop  = 171.;
  //double pTmin = 100.;

  // Generator. Shorthand for the event and particle data.
  Pythia pythia;
  Event& event = pythia.event;
  ParticleData& pData = pythia.particleData;

  // Set up beams.
  pythia.settings.parm("Beams:eCM", eCM);
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("PDF:lepton = off");

  // Set up top pair production. Require pTmin kinematics cut.
  // pythia.readString("Top:gg2ttbar = on");
  // pythia.readString("Top:qqbar2ttbar = on");
  // pythia.settings.parm("PhaseSpace:pTHatMin", pTmin);
  pythia.readString("Top:ffbar2ttbar(s:gmZ) = on");

  // Set top mass and width. No resonance description.
  pData.m0( 6, mTop);
  pData.doForceWidth( 6, true);
  pData.mWidth( 6, 2.);

  // Set up user hook.
  MyUserHooks* myUserHooks = new MyUserHooks();
  pythia.setUserHooksPtr(myUserHooks);

  // Switch off key components. Useful for first checks, but not full run.
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("HadronLevel:Hadronize = off");

  // Possibility to switch off particle data and event listings.
  //pythia.readString("Init:showChangedSettings = off");
  //pythia.readString("Init:showChangedParticleData = off");
  //pythia.readString("Next:numberShowInfo = 0");
  //pythia.readString("Next:numberShowProcess = 0");
  //pythia.readString("Next:numberShowEvent = 0");

  // Useful options: alternative tunes, modifed top fragmentation.
  //pythia.readString("Tune:ee = 7");
  //pythia.readString("Tune:pp = 14");
  //pythia.readString("StringZ:rFactH = 1.");

  // Debug options.
  pythia.readString("Check:nErrList = 5");

  // Initialize.
  cout << " top width before init = " << pData.mWidth(6) << endl;
  pythia.init();
  cout << " top width after  init = " << pData.mWidth(6) << endl;

  // Histograms.
  Hist nChargedH("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate events. Skip. if failure.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Loop over final charged particles in the event.
    int nCharged = 0;
    for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
    nChargedH.fill( nCharged );

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();
  cout << nChargedH << egBef << egAft;

  // Done.
  return 0;
}
