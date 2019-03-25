#ifndef Low_Energy_Resonance_H
#define Low_Energy_Resinance_H

#include "Pythia8/Event.h"
#include "Pythia8/LowEnergyData.h"

namespace Pythia8 {

class LowEnergyResonance {
public:

  void initPtr(Rndm* rndmPtrIn, ParticleData* particleDataPtrIn, LowEnergyData* lowEnergyDataPtrIn)
  { rndmPtr = rndmPtrIn; particleDataPtr = particleDataPtrIn; lowEnergyDataPtr = lowEnergyDataPtrIn; }

  bool collide(int i1, int i2, Event& event);

  
  double getPartialResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;

  double getPartialElasticResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;

  double getResonanceSigma(int idA, int idB, double eCM) const;

  double getElasticResonanceSigma(int idA, int idB, double eCM) const;


private:

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  LowEnergyData* lowEnergyDataPtr;

};

}

#endif