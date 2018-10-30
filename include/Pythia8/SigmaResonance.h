#ifndef Pythia8_SigmaResonance_H
#define Pythia8_SigmaResonance_H

#include "Pythia8/Interpolator.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/ResonanceData.h"

namespace Pythia8 {

class SigmaResonance {

public:

  SigmaResonance() { };

  void initPtr(Rndm* rndmPtrIn, ParticleData* particleDataPtrIn, 
    ResonanceData* resonanceDataPtrIn) {
    rndmPtr = rndmPtrIn; particleDataPtr = particleDataPtrIn;
    resonanceDataPtr = resonanceDataPtrIn;
  }

  double sigma(int idA, int idB, double eCM) const;

  pair<int, int> pickProducts(int idAIn, int idBIn, double eCM);

  const Interpolator& sigmaDistribution(int idA, int idB) const; 
  const Interpolator& sigmaDistribution(int idA, int idB, int idC, int idD) const;

private:

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  ResonanceData* resonanceDataPtr;


};

}

#endif