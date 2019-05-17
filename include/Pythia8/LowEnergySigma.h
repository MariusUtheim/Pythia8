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

  void init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, ParticleWidths* particleWidthsPtrIn);


  // Get the total cross section for the specified collision
  double sigmaTotal(int idA, int idB, double eCM) const;

  // Get the partial cross section for the specified collision and process.
  // @TODO: Decide on exact scheme for process ids
  // 1: non-diffractive | 2: elastic | 3: SD (XB) | 4: SD (AX)
  // 5: DD | 6: annihilation | 7: resonant | 8: excitation
  // >100: resonant through the particle specified by the id
  double sigmaPartial(int idA, int idB, double eCM, int process) const;

  // Gets all partial cross sections for the specified collision. 
  // This is used when all cross sections are needed to determine which 
  // process to execute. The keys are process ids, the values are cross 
  // sections, which give the relative frequency of that process.
  map<int, double> sigmaPartial(int idA, int idB, double eCM) const;

private:

  Info* infoPtr;

  ParticleData* particleDataPtr;

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
  double sigmaTotalBB(int idA, int idB, double eCM) const;
  double sigmaElasticBB(int idA, int idB, double eCM) const;
  double sigmaStrEx(int idA, int idB, double eCM) const;
  map<int, double> sigmaPartialBB(int idA, int idB, double eCM) const;


  // @TODO: Rethink the interface here

  // BBbar
  double sigmaTotalBBbar(int idA, int idB, double eCM) const;
  map<int, double> sigmaPartialBBbar(int idA, int idB, double eCM) const;
  
  // XM
  ParticleWidths* particleWidthsPtr;
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