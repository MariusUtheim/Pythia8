#ifndef Pythia8_LowEnergySigma_H
#define Pythia8_LowEnergySigma_H

#include "Pythia8/Info.h"
#include "Pythia8/LowEnergyResonance.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

class LowEnergySigma {
public:

  void initPtr(Info* infoPtrIn, ParticleData* particleDataPtrIn,
    LowEnergyResonance* lowEnergyResPtrIn) {
    infoPtr = infoPtrIn; particleDataPtr = particleDataPtrIn; 
    lowEnergyResPtr = lowEnergyResPtrIn;
  }

  double sigmaTotal(int idA, int idB, double eCM) const;

  // Scattering only through the specified process
  // 1: non-diffractive | 2: elastic | 3: SD (XB) | 4: SD (AX)
  // 5: DD | 6: annihilation | 7: resonant
  // >100: resonant through the particle specified by the id
  double sigmaPartial(int idA, int idB, double eCM, int process) const;

  // Gets all possible cross sections for the specified particles. Entry 0
  // gives total cross section. Depending on the collision class, all the 
  // following keys will be in the resulting map:
  //   BB:    0, 1, 2, 3, 4, 5
  //   BBbar: 0, 1, 2, 3, 4, 5, 6
  //   XM:    0, 1, 2, 7, plus each possible resonance id
  map<int, double> sigmaPartial(int idA, int idB, double eCM) const;

private:

  Info* infoPtr;

  ParticleData* particleDataPtr;

  LowEnergyResonance* lowEnergyResPtr;

  double aqm(int idA, int idB) const;
  double aqmNN() const;

  double sigmaTotalBB(int idA, int idB, double eCM) const;

  double sigmaElasticBB(int idA, int idB, double eCM) const;

  double sigmaTotalBBbar(int idA, int idB, double eCM) const;
  map<int, double> sigmaPartialBBbar(int idA, int idB, double eCM) const;
  
  double sigmaTotalXM(int idX, int idM, double eCM) const;
  double sigmaElasticXM(int idX, int idM, double eCM) const;
  double sigmaNondiffXM(int idX, int idM, double eCM) const;

  bool hasParametrisation(int idA, int idB) const;

  // Select process type:
  //  1: baryon-baryon or antibaryon-antibaryon
  //  2: baryon-antibaryon
  //  3: hadron-meson
  int selectProcess(int idA, int idB) const;

  int canonicalForm(int& idA, int& idB) const;
  int strangenessFactor(int id) const {
    return 1. - 0.4 * particleDataPtr->strangeness(abs(id))
                    / (particleDataPtr->isBaryon(id) ? 3. : 2.);
  }
};

}

#endif