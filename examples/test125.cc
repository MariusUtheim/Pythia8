// test125.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Simple illustration how to implement Dark Matter pair production at the LHC.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// A derived class for q qbar -> DMmediator -> DM DM.

class Sigma1qqbar2DMmediator : public Sigma1Process {

public:

  // Constructor.
  Sigma1qqbar2DMmediator() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). Assumed flavour-independent so simple.
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles. Not needed for isotropic decays.
  // virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "q qbar -> DMmediator";}
  virtual int    code()       const {return 10001;}
  virtual string inFlux()     const {return "qqbarSame";}
  virtual int    resonanceA() const {return idDMmed;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    idDMmed;
  double mRes, GammaRes, m2Res, GamMRat, normDMmed2qqbar, sigma;

  // Pointer to properties of DMmediator, to access decay width.
  ParticleDataEntry* particlePtr;

};

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1qqbar2DMmediator::initProc() {

  // Store DMmediator mass and width for propagator.
  idDMmed  = 54;
  mRes     = particleDataPtr->m0(idDMmed);
  GammaRes = particleDataPtr->mWidth(idDMmed);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Arbitrary coupling strength normalization; inn principle related to width.
  normDMmed2qqbar = 0.0001;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(idDMmed);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat); first step when inflavours unknown.

void Sigma1qqbar2DMmediator::sigmaKin() {

  // Incoming width with colour factor.
  double widthIn  = normDMmed2qqbar * mH / 3.;

  // Breit-Wigner, including some (guessed) spin factors.
  double sigBW    = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Outgoing width: only includes channels left open.
  double widthOut = particlePtr->resWidthOpen(idDMmed, mH);

  // Total answer.
  sigma = widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1qqbar2DMmediator::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idDMmed);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 1, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

int main() {

  // Number of events to generate. Max number of errors.
  int nEvent = 1000;
  int nAbort = 5;

  // Pythia generator.
  Pythia pythia;

  // Set DM and DMmediator mass and width.
  pythia.readString("51:m0 = 112.");
  pythia.readString("54:m0 = 222.");
  pythia.readString("54:mWidth = 2.");

  // Allow DMmediator to decay to DM(s=0) DM(s=0) pair only.
  // Note: do not want to use 54:onIfMatch = 51 51 since onIfMatch
  // is sign-insensitive so also switches on 51 -51.
  pythia.readString("54:onMode = off");
  pythia.readString("54:19:onMode = on");

  // Create instance of a class to generate the q qbar -> DMmediator process
  // from an external matrix element. Hand in pointer to Pythia.
  SigmaProcess* sigma1DMmediator = new Sigma1qqbar2DMmediator();
  pythia.setSigmaPtr(sigma1DMmediator);

  // Optionally only compare cross sections.
  //pythia.readString("PartonLevel:all = off");
  pythia.readString("Check:nErrList = 2");

  // Initialization for LHC.
  pythia.init();

  // Book histogram.
  Hist mDMmed("DMmediator mass", 100, 0., 400.);
  Hist pTDMmed("DMmediator pT",  100, 0., 400.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Find last DMmediator copy.
    int iDMmed = 0;
    for (int i = 1; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 54) iDMmed = i;
    if (iDMmed == 0) continue;

    // DMmediator mass and pT. End of event loop.
    mDMmed.fill( pythia.event[iDMmed].m() );
    pTDMmed.fill( pythia.event[iDMmed].pT() );
  }

  // Final statistics. Print histograms.
  pythia.stat();
  cout << mDMmed << pTDMmed;

  // Done.
  delete sigma1DMmediator;
  return 0;
}
