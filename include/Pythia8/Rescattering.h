#ifndef Pythia8_Rescattering_H
#define Pythia8_Rescattering_H

#include "Pythia8/Basics.h"
#include "Pythia8/CrossSectionData.h"
#include "Pythia8/Event.h"
#include "Pythia8/ParticleDecays.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/ResonanceDecays.h"

namespace Pythia8 {

class Rescattering {

public:

  Rescattering() {}

  void init(Info* infoPtrIn, Rndm* rndmPtrIn)
  { infoPtr = infoPtrIn; rndmPtr = rndmPtrIn; }

  // @TODO calculate origin in rescatter call?
  void rescatter(int idA, int idB, Vec4 origin, Event& event);

  bool calcRescatterOrigin(int idA, int idB, Event& event, Vec4& originOut);

private: 

  Info* infoPtr;

  Rndm* rndmPtr;

};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
