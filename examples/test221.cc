// main62.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how you can use UserHooks to set angular decay distributions
// for undecayed resonances from Les Houches input using the polarization
// information of the boson defined in its rest frame.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================
// Book a histogram to test the angular distribution in the UserHook.
// It is booked here so that it is global.

Hist cosRaw("cos raw",100,-1.,1.);

//==========================================================================

// Write own derived UserHooks class.
// Assumptions in this particular case:
// The W+- bosons were undecayed in the Les Houches Events input file,
// and subsequently decayed isotropically by the Pythia machinery.
// Now the angular distribution will be corrected for each W,
// based on the polarization value stored in the LHEF.
// For W- this is (1 -+ cos(theta))^2 for +-1, sin^2(theta) for 0,
// and isotropic for 9. For W+ it is flipped (i.e. theta->pi-theta).
// The Pythia decay products (i.e. the branching ratios) are retained.

class MyUserHooks : public UserHooks {

public:

  // Constructor and destructor do nothing.
  MyUserHooks() {}
  ~MyUserHooks() {}

  // Allow a veto for the process level, to gain access to decays.
  bool canVetoProcessLevel() {return true;}

  // Access the event after resonance decays.
  bool doVetoProcessLevel(Event& process) {

    // Identify decayed W+- bosons for study.
    // Assume isotropic decay if polarization is unphysically big (|pol|>2)
    for (int i = 0; i < process.size(); ++i) {
      if (process[i].idAbs() == 23 && process[i].status() == -22) {

	// Identify W+- daughters, properly ordered.
	int idV = process[i].id();
	int i1  = process[i].daughter1();
	int i2  = process[i].daughter2();

	double theta = acos( 2. * rndmPtr->flat() - 1.);
	double cosThe = cos(theta);
	double sinThe = sin(theta);
	double phi   = 2. * M_PI * rndmPtr->flat();

	process[i1].id(443);
	process[i2].id(22);
	process[i1].col(0);
	process[i1].acol(0);
	process[i2].col(0);
	process[i2].acol(0);
	process[i1].status(82);
	process[i2].status(1);
	// Set up decay in rest frame of W+-.
	double mV = process[i].m();
	double m1 = particleDataPtr->m0(443);
	double m2 = 0.0;
	// energy of first decay product in W rest frame
	double e1 = 0.5* (pow2(mV) + pow2(m1) - pow2(m2))/mV;
	// momentum magnitude of the same particle
	double pA = sqrt(pow2(e1) - pow2(m1));
	Vec4 p1( pA * sinThe * cos(phi), pA *sinThe * sin(phi), pA * cosThe, e1);
	Vec4 p2   = Vec4(0,0,0,mV) - p1;

	// Angular reference axis defined as opposite the mother
	// direction in W rest frame.
	Vec4 pM = process[process[i].mother1()].p();
	pM.bstback( process[i].p() );
	pM.flip3();
	RotBstMatrix Mrotbst;
	Mrotbst.rot( pM);

	// Rotate and boost W decay products.
	Mrotbst.bst( process[i].p() );
	p1.rotbst(Mrotbst);
	p2.rotbst(Mrotbst);
	process[i1].p( p1 );
	process[i2].p( p2 );
	process[i1].m( p1.mCalc() );
	process[i2].m( p2.mCalc() );
      }
      // End of loop over W's. Do not veto any events.
    }
    return false;
  }

  // Select polar angle for the W decay.
  double selectAngle( double inputSpin, double inputId ) {

    // Set up initial angles.
    double rdNow = rndmPtr->flat();
    double cosThe;
    // Small number to distinguish -1,1, and 0 with round-off
    double eps = 1e-10;

    // W+ distribution is "opposite" of W-.
    if (inputId > 0) inputSpin *= -1;

    // Different decay angular distributions.
    // 3/8 * (1 - cos(theta))^2  ++
    if (inputSpin > eps) {
      cosThe = max( 1.0 - 2.0 * pow(rdNow, 1./3.), -1.0);
    // 3/8 * (1 + cos(theta))^2  --
    } else if (inputSpin < -eps) {
      cosThe = min( 2.0 * pow(rdNow, 1./3.) - 1.0,  1.0);
    // 3/4 * sin(theta)^2        00
    // Solution of cubic equation that yields the correct result.
    } else {
      double theA = (acos(1.0 - 2.0 * rdNow) + 4.0 * M_PI) / 3.0;
      cosThe = 2.0 * cos(theA);
    }

    // Return the selected cos(theta) value.
    return cosThe;
  }

};

//==========================================================================

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up to do a user veto and send it in. Initialize.
  MyUserHooks* myUserHooks = new MyUserHooks();
  pythia.setUserHooksPtr( myUserHooks);
  pythia.readFile("test221.cmnd");
  pythia.init();

  // Histograms.
  Hist cosPlus(" cos(the) W- -> f",100,-1.0,1.0);
  Hist cosMinus(" cos(the) W+ -> fbar",100,-1.0,1.0);
  Hist polarization(" polarization",10,-2.0,2.0);
  Hist energy(" energy ",100,0.0,100.0);

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  //  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events.
    pythia.next();

    // Loop through event, looking for a W and its daughters.
    for (int i = 0; i < event.size(); ++i) {
      // Select W boson when it decays to two partons, not when it is
      // a recoil in FSR.
      if (event[i].idAbs() == 24
	  && event[i].daughter1() != event[i].daughter2() ) {
	int i1 = event[i].daughter1();
	// Angular distribution is defined with respect to the decay product
	// with the same sign charge as the W boson
	if (event[i1].id() * event[i].id() > 0 ) i1 = event[i].daughter2();

	// Reconstruct W+- decay angle. Histogram information.
	Vec4 p1 = event[i1].p();
	RotBstMatrix Mrotbst;
	Mrotbst.bst( event[i].p(), Vec4( 0., 0., 0., event[i].m()) );
	p1.rotbst( Mrotbst );
	energy.fill( p1.e() );
	// Boost mother to rest frame of W and flip direction
	Vec4 p2 = event[event[i].mother1()].p();
	p2.rotbst( Mrotbst );
	p2.flip3();
	double costhe = costheta( p1, p2 );
	polarization.fill( event[i].pol() );
	if( event[i].id() > 0 ) cosPlus.fill( costhe );
	else                   cosMinus.fill( costhe );
      }
    }

  // End of event loop.
  }

  // Statistics. Histograms.
  pythia.stat();
  cout << cosPlus << cosMinus << polarization << energy << cosRaw;

  // Done.
  delete myUserHooks;
  return 0;
}
