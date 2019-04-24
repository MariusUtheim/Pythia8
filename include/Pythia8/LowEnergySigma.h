#ifndef Pythia8_LowEnergySigma_H
#define Pythia8_LowEnergySigma_H

#include "Pythia8/Info.h"
#include "Pythia8/LowEnergyResonance.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

//==========================================================================

// Gets cross sections for all hadron-hadron collisions at low energies

class LowEnergySigma {
public:

  void initPtr(Info* infoPtrIn, ParticleData* particleDataPtrIn,
    LowEnergyResonance* lowEnergyResPtrIn) {
    infoPtr = infoPtrIn; particleDataPtr = particleDataPtrIn; 
    lowEnergyResPtr = lowEnergyResPtrIn;
  }

  // Get the total cross section for the specified collision
  double sigmaTotal(int idA, int idB, double eCM) const;

  // Get the partial cross section for the specified collision and process.
  // 1: non-diffractive | 2: elastic | 3: SD (XB) | 4: SD (AX)
  // 5: DD | 6: annihilation | 7: resonant
  // >100: resonant through the particle specified by the id
  double sigmaPartial(int idA, int idB, double eCM, int process) const;

  // Gets all partial cross sections for the specified collision. 
  // This is used when all cross sections are needed to determine which 
  // process to execute.
  // Entry 0 gives total cross section. Depending on the collision class, all 
  // the following keys will be in the resulting map:
  //   BB:    0, 1, 2, 3, 4, 5
  //   BBbar: 0, 1, 2, 3, 4, 5, 6
  //   XM:    0, 1, 2, 7, plus each possible resonance id
  map<int, double> sigmaPartial(int idA, int idB, double eCM) const;

private:

  Info* infoPtr;

  ParticleData* particleDataPtr;

  LowEnergyResonance* lowEnergyResPtr;

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

  // @TODO: Rethink the interface here

  // BBbar
  double sigmaTotalBBbar(int idA, int idB, double eCM) const;
  map<int, double> sigmaPartialBBbar(int idA, int idB, double eCM) const;
  
  // XM
  double sigmaTotalXM(int idX, int idM, double eCM) const;
  double sigmaElasticXM(int idX, int idM, double eCM) const;
  double sigmaInelXM(int idX, int idM, double eCM) const;

  // Select process type:
  //  1: baryon-baryon or antibaryon-antibaryon
  //  2: baryon-antibaryon
  //  3: hadron-meson
  int selectProcess(int idA, int idB) const;


};

}

#endif