// test164.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class.

class MyUserHooks : public UserHooks {

public:

  // Constructor and destructor.
  MyUserHooks() {}
  ~MyUserHooks() {}

  // Allow process level event to be modified.
  virtual bool canVetoProcessLevel() {return true;}

  // ...which gives access to the event after process level has been set.
  virtual bool doVetoProcessLevel(Event& process) {

    // Modify flavours, colours and masses.
    process[6].id(421);
    process[7].id(22);
    process[6].cols( 0, 0);
    process[7].cols( 0, 0);
    double mD   = particleDataPtr->mSel(421);
    process[6].m( mD);
    process[7].m( 0.);

    // Modify kinematics: isotropic decay in rest frame and then boost.
    double mZ   = process[5].m();
    double eD   = (mZ * mZ + mD * mD) / (2. * mZ);
    double pD   = (mZ * mZ - mD * mD) / (2. * mZ);
    double cThe = 2. * rndmPtr->flat() - 1.;
    double sThe = sqrt(1. - cThe * cThe);
    double phi  = 2. * M_PI * rndmPtr->flat();
    double cPhi = cos(phi);
    double sPhi = sin(phi);
    process[6].p(  pD * sThe * cPhi,  pD * sThe * sPhi,  pD * cThe, eD);
    process[7].p( -pD * sThe * cPhi, -pD * sThe * sPhi, -pD * cThe, pD);
    process[6].bst( process[5].p() );
    process[7].bst( process[5].p() );

    // Return with modified event and no veto.
    return false;
  }

};

//==========================================================================

int main() {

  // Number of events.
  int nEvent = 1000;

  // Generator and process.
  Pythia pythia;
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

  // Set up to do a user veto and send it in.
  MyUserHooks* myUserHooks = new MyUserHooks();
  pythia.setUserHooksPtr( myUserHooks);

  // Set D0 decay.
  pythia.readString("421:onMode = off");
  pythia.readString("421:onIfMatch = 321 211");

  // Initialize. Histogram.
  pythia.init();
  Hist mZ("Z0 mass", 100, 0., 200.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    mZ.fill( pythia.process[5].m() );

  // End of event loop. Statistics and histogram.
  }
  pythia.stat();
  cout << mZ;

  // Done.
  return 0;
}
