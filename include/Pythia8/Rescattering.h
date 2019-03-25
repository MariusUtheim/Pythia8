#ifndef Pythia8_Rescattering_H
#define Pythia8_Rescattering_H

#include "Pythia8/Basics.h"
#include "Pythia8/LowEnergyHadHad.h"
#include "Pythia8/Event.h"

namespace Pythia8 {

class Rescattering {

public:

  Rescattering() {}

  void init(Info* infoPtrIn, Rndm* rndmPtrIn, ParticleData* particleDataPtrIn,
    LowEnergyHadHad* resDataPtrIn)
  { infoPtr = infoPtrIn; rndmPtr = rndmPtrIn; 
  particleDataPtr = particleDataPtrIn; resDataPtr = resDataPtrIn; }

  // @TODO calculate origin in rescatter call?
  void rescatter(int idA, int idB, Vec4 origin, Event& event);

  bool calcRescatterOrigin(int idA, int idB, Event& event, Vec4& originOut);

private: 

  Info* infoPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  // @TODO Better name
  LowEnergyHadHad* resDataPtr;
};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
