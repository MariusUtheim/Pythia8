#ifndef Pythia8_LowEnergySigma_H
#define Pythia8_LowEnergySigma_H

#include "Pythia8/Info.h"
#include "Pythia8/ParticleWidths.h"
#include "Pythia8/ParticleData.h"
#include "SigmaTotal.h"

namespace Pythia8 {

//==========================================================================

// Gets cross sections for all hadron-hadron collisions at low energies

class LowEnergySigma {
public:

  void init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn, ParticleWidths* particleWidthsPtrIn);


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
  //map<int, double> sigmaPartial(int idA, int idB, double eCM) const; 
  
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

  mutable SigmaSaSDL sigmaSaSDL;

  // Orders the two inputs in a canonical way, and returns an id for the 
  // collision type (see .cc file for a detailed specification)
  // 1: BB |  2: BBbar |  3: XM
  int canonicalForm(int& idA, int& idB) const;

  // Get cross section according to additive quark model.
  // Used for generic BB collisions, and for scale factors
  double aqm(int idA, int idB) const;
  // Get AQM nucleon-nucleon cross section (always 40 mb). For scale factors
  double aqmNN() const;

  // BB
  double BBTotal(int idA, int idB, double eCM) const;
  double BBElastic(int idA, int idB, double eCM) const;
  double BBStrEx(int idA, int idB, double eCM) const;
  double BBNonDiff(int idA, int idB, double eCM) const;
  double BBDiffractiveAX(int idA, int idB, double eCM) const;
  double BBDiffractiveXB(int idA, int idB, double eCM) const;
  double BBDiffractiveXX(int idA, int idB, double eCM) const;
  double BBExcite(int idA, int idB, double eCM) const;

  map<int, double> sigmaPartialBB(int idA, int idB, double eCM) const;


  // @TODO: Rethink the interface here

  // BBbar
  double BBbarTotal(int idA, int idB, double eCM) const;
  double BBbarElastic(int idA, int idB, double eCM) const;
  double BBbarNonDiff(int idA, int idB, double eCM) const;
  double BBbarDiffractiveAX(int idA, int idB, double eCM) const;
  double BBbarDiffractiveXB(int idA, int idB, double eCM) const;
  double BBbarDiffractiveXX(int idA, int idB, double eCM) const;
  double BBbarAnnihilation(int idA, int idB, double eCM) const;

  map<int, double> sigmaPartialBBbar(int idA, int idB, double eCM) const;
  
  // XM
  ParticleWidths* particleWidthsPtr;
  double XMTotal(int idX, int idM, double eCM) const;
  double XMElastic(int idX, int idM, double eCM) const;
  double XMResonant(int idX, int idM, double eCM) const;


  double sigmaTotalXM(int idX, int idM, double eCM) const;
  double sigmaElasticXM(int idX, int idM, double eCM) const;
  double sigmaResTotalXM(int idX, int idM, double eCM) const;
  double sigmaResPartialXM(int idX, int idM, int idR, double eCM) const;
  map<int, double> sigmaResXM(int idX, int idM, double eCM) const;
  double sigmaStringXM(int idX, int idM, double eCM) const;
  // @TODO: Make a more intutive system?
  // The signature of a particle is the three digit number BQS, where B is 
  // baryon number, Q is charge signature and S is strangeness signature.
  // A resonance can be formed only if it conserves the total signature. 
  // The charge signature of a particle with charge q is given by chargeType if
  // charge is positive and 10 + chargeType if it is negative. This way, charge
  // signature is always positive. Strangeness signature is defined similarly.
  map<int, vector<int>> signatureToParticles;

};

}

#endif