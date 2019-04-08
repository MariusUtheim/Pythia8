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
  // >101: resonant through the particle specified by the id
  double sigmaPartial(int idA, int idB, double eCM, int process) const;

private:

  Info* infoPtr;

  ParticleData* particleDataPtr;

  LowEnergyResonance* lowEnergyResPtr;

  double sigmaTotalBB(int idA, int idB, double eCM) const;

  double sigmaTotalBBbar(int idA, int idB, double eCM) const;

  double sigmaTotalXM(int idX, int idM, double eCM) const;

  bool hasParametrisation(int idA, int idB) const;

  // Select process type:
  //  1: baryon-baryon or antibaryon-antibaryon
  //  2: baryon-antibaryon
  //  3: hadron-meson
  int selectProcess(int idA, int idB) const;

  int strangenessFactor(int id) const {
    return particleDataPtr->strangeness(abs(id))
         / (particleDataPtr->isBaryon(id) ? 3. : 2.);
  }
};

}

#endif