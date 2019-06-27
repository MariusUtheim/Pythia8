#ifndef Pythia8_LowEnergySigma_H
#define Pythia8_LowEnergySigma_H

#include "Pythia8/Info.h"
#include "Pythia8/HadronWidths.h"
#include "Pythia8/ParticleData.h"
#include "SigmaTotal.h"

namespace Pythia8 {

//==========================================================================

// Gets cross sections for all hadron-hadron collisions at low energies

class LowEnergySigma {
public:

  void init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn, HadronWidths* hadronWidthsPtrIn);


  // Get the total cross section for the specified collision
  double sigmaTotal(int idA, int idB, double eCM) const;

  // Get the partial cross section for the specified collision and process.
  // type | 0: mix; | 1: nondiff; | 2 : el; | 3: SD (XB); | 4: SD (AX);
  //      | 5: DD;  | 6: CD (AXB, not implemented) 
  //      | 7: excitation | 8: annihilation | 9: resonant
  //      | >100: resonant through the specified resonance particle
  double sigmaPartial(int idA, int idB, double eCM, int type) const;

  // Gets all partial cross sections for the specified collision. 
  // This is used when all cross sections are needed to determine which 
  // process to execute. Returns false if no processes are available.
  bool sigmaPartial(int idA, int idB, double eCM, 
    vector<int>& typesOut, vector<double>& sigmasOut) const;

  // Picks a process randomly according to their partial cross sections
  int pickProcess(int idA, int idB, double eCM);

  // Picks a resonance according to their partial cross sections
  int pickResonance(int idA, int idB, double eCM);

private:

  Info* infoPtr;

  ParticleData* particleDataPtr;

  Rndm* rndmPtr;

  // @TODO Just put the actual formulas in there
  mutable SigmaSaSDL sigmaSaSDL;

  // Orders the two inputs in a canonical way, and returns an id for the 
  // collision type (see .cc file for a detailed specification)
  // 1: BB |  2: BBbar |  3: XM |  4: -XM
  int canonicalForm(int& idA, int& idB) const;

  // Get cross section according to additive quark model.
  // Used for generic BB collisions, and for scale factors
  double aqm(int idA, int idB) const;
  // Get AQM nucleon-nucleon cross section (always 40 mb). For scale factors
  double aqmNN() const;

  // BB partial cross sections
  double BBTotal(int idA, int idB, double eCM) const;
  double BBElastic(int idA, int idB, double eCM) const;
  double BBStrEx(int idA, int idB, double eCM) const;
  double BBNonDiff(int idA, int idB, double eCM) const;
  double BBDiffractiveAX(int idA, int idB, double eCM) const;
  double BBDiffractiveXB(int idA, int idB, double eCM) const;
  double BBDiffractiveXX(int idA, int idB, double eCM) const;
  double BBExcite(int idA, int idB, double eCM) const;

  // BBbar partial cross sections
  double BBbarTotal(int idA, int idB, double eCM) const;
  double BBbarElastic(int idA, int idB, double eCM) const;
  double BBbarNonDiff(int idA, int idB, double eCM) const;
  double BBbarDiffractiveAX(int idA, int idB, double eCM) const;
  double BBbarDiffractiveXB(int idA, int idB, double eCM) const;
  double BBbarDiffractiveXX(int idA, int idB, double eCM) const;
  double BBbarAnnihilation(int idA, int idB, double eCM) const;

  // XM partial cross sections
  HadronWidths* hadronWidthsPtr;
  double XMTotal(int idX, int idM, double eCM) const;
  double XMNonDiffractive(int idX, int idM, double eCM) const;
  double XMElastic(int idX, int idM, double eCM) const;
  double XMResonant(int idX, int idM, double eCM) const;
  double XMResonantPartial(int idX, int idM, int idR, double eCM) const;
  vector<int> possibleResonances(int idX, int idM) const;

  // @TODO: Make a more intutive system?
  // The signature of a particle is the three digit number BQS, where B is 
  // baryon number, Q is charge signature and S is number of strange quarks.
  // A resonance can be formed only if it conserves the total signature. 
  // The charge signature of a particle with charge q is given by chargeType if
  // charge is positive and 10 + chargeType if it is negative. This way, charge
  // signature is always positive. Strangeness signature is defined similarly.
  map<int, vector<int>> signatureToParticles;
  int getSignature(int baryon, int charge, int nStrange) const;

};

}

#endif