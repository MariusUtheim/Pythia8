#ifndef Low_Energy_Diffractive_H
#define Low_Energy_Diffractive_H

#include "Pythia8/Event.h"
#include "Pythia8/LowEnergyData.h"
#include "Pythia8/LowEnergyProcess.h"

namespace Pythia8 {

class LowEnergyDiffractive : public LowEnergyProcess {
public:

  void initPtr(Rndm* rndmPtrIn, ParticleData* particleDataPtrIn, LowEnergyData* lowEnergyDataPtrIn)
  { rndmPtr = rndmPtrIn; particleDataPtr = particleDataPtrIn; lowEnergyDataPtr = lowEnergyDataPtrIn; }

  bool collide(int i1, int i2, Event& event) const;

  double getDiffractiveSigma(int idA, int idB, double eCM) const;


private:

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  LowEnergyData* lowEnergyDataPtr;
};

}

#endif