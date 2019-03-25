// test158.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program, loosely based on main17.
// It illustrates external sequential decays D+ -> K*0bar pi+ -> K- pi+ pi+.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// A derived class to do D+ decays.

class DplusDecay : public DecayHandler {

public:

  // Constructor.
  DplusDecay(ParticleData* pdtPtrIn, Rndm* rndmPtrIn) {times = 0;
    pdtPtr = pdtPtrIn; rndmPtr = rndmPtrIn;}

  // Routine for doing the sequantial decay.
  bool chainDecay(vector<int>& idProd, vector<int>& motherProd,
    vector<double>& mProd, vector<Vec4>& pProd, int iDec, const Event& event);

private:

  // Count number of times DplusDecay is called.
  int times;

  // Pointer to the particle data table.
  ParticleData* pdtPtr;

  // Pointer to the random number generator.
  Rndm* rndmPtr;

};

//--------------------------------------------------------------------------

// The actual D+ decay routine.
// Not intended for realism, just to illustrate the principles.

bool DplusDecay::chainDecay(vector<int>& idProd,  vector<int>& motherProd,
  vector<double>& mProd, vector<Vec4>& pProd, int iDec, const Event& ) {

  // Do decay D+ -> K*0bar pi+; store the products.
  int iSgn = (idProd[0] > 0) ? 1 : -1;
  idProd.push_back(-313 * iSgn);
  idProd.push_back(211 * iSgn);
  motherProd.push_back(0);
  motherProd.push_back(0);

  // Decay masses, here from Pythia tables, also stored.
  mProd.push_back(pdtPtr->m0(313));
  mProd.push_back(pdtPtr->m0(211));

  // Calculate muon energy and momentum in D+ rest frame.
  double m0s = pow2(mProd[0]);
  double m1s = pow2(mProd[1]);
  double m2s = pow2(mProd[2]);
  double pAbs = 0.5 * sqrt( pow2(m0s - m1s - m2s) - 4. * m1s * m2s) / mProd[0];

  // Assume decay angles isotropic in rest frame.
  double cosTheta = 2. * rndmPtr->flat() - 1.;
  double sinTheta = sqrt(max(0., 1. - cosTheta * cosTheta));
  double phi = 2. * M_PI * rndmPtr->flat();
  double px = pAbs * sinTheta * cos(phi);
  double py = pAbs * sinTheta * sin(phi);
  double pz = pAbs * cosTheta;

  // Define four-vectors in the D+ rest frame.
  Vec4 p1(  px,  py,  pz, 0.5 * (m0s + m1s - m2s) / mProd[0]);
  Vec4 p2( -px, -py, -pz, 0.5 * (m0s + m2s - m1s) / mProd[0]);

  // Boost them by velocity vector of the D+ mother and store.
  p1.bst(pProd[0]);
  p2.bst(pProd[0]);
  pProd.push_back(p1);
  pProd.push_back(p2);

  // Do decay K*0bar -> K- pi+; store the products.
  idProd.push_back(-321 * iSgn);
  idProd.push_back(211 * iSgn);
  motherProd.push_back(1);
  motherProd.push_back(1);

  // Decay masses, here from Pythia tables, also stored.
  mProd.push_back(pdtPtr->m0(321));
  mProd.push_back(pdtPtr->m0(211));

  // Calculate muon energy and momentum in D+ rest frame.
  double m3s = pow2(mProd[3]);
  double m4s = pow2(mProd[4]);
  pAbs = 0.5 * sqrt( pow2(m1s - m3s - m4s) - 4. * m3s * m4s) / mProd[1];

  // Assume decay angles isotropic in rest frame.
  cosTheta = 2. * rndmPtr->flat() - 1.;
  sinTheta = sqrt(max(0., 1. - cosTheta * cosTheta));
  phi = 2. * M_PI * rndmPtr->flat();
  px = pAbs * sinTheta * cos(phi);
  py = pAbs * sinTheta * sin(phi);
  pz = pAbs * cosTheta;

  // Define four-vectors in the K*0bar rest frame.
  Vec4 p3(  px,  py,  pz, 0.5 * (m1s + m3s - m4s) / mProd[1]);
  Vec4 p4( -px, -py, -pz, 0.5 * (m1s + m4s - m3s) / mProd[1]);

  // Boost them by velocity vector of the K*0bar mother and store.
  p3.bst(pProd[1]);
  p4.bst(pProd[1]);
  pProd.push_back(p3);
  pProd.push_back(p4);

  // Print message the first few times, to show that it works.
  if (times++ < 10) cout << "\n External D+ decay performed, D+ in line "
    << iDec << "\n";

  // Done
  return true;

}

//==========================================================================

int main() {

  // Number of events to generate and to list. Max number of errors.
  int nEvent = 10000;
  int nList  = 2;
  int nAbort = 5;

  // Pythia generator.
  Pythia pythia;

  // Initialization for charm production. Reduce event size.
  pythia.readString("HardQCD:qqbar2ccbar = on");
  pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");

  // Set fictitious long K* lifetime to check secondary vertices.
  pythia.readString("313:tau0 = 5.");

  // A class to do D+ decays externally.
  DecayHandler* handleDecays = new DplusDecay(&pythia.particleData,
    &pythia.rndm);

  // The list of particles the class can handle.
  vector<int> handledParticles;
  handledParticles.push_back(411);

  // Hand pointer and list to Pythia.
  pythia.setDecayPtr( handleDecays, handledParticles);

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialization.
  pythia.init();

  // Book histograms.
  Hist pThard("pTHat of hard subprocess", 100, 0., 50.);
  Hist pTJPsi("pT of J/Psi", 100, 0., 50.);

  // Begin event loop.
  int iList = 0;
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Histogram pThard spectrum of process.
    double pTHat = pythia.info.pTHat();
    pThard.fill( pTHat );

    // Look for event with externally handled decays.
    bool externalDecay = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      int status = pythia.event[i].statusAbs();
      if (status == 93 || status == 94) {externalDecay = true; break;}
    }

    // List first few events with external decay.
    if (externalDecay && ++iList <= nList) pythia.event.list(true);

  // End of event loop.
  }

  // Final statistics. Print histograms.
  pythia.stat();
  cout << pThard;

  // Done.
  delete handleDecays;
  return 0;
}
