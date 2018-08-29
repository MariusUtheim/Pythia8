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

  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn, CrossSectionData* crossSectionDataPtrIn,
    UserHooks* userHooksIn) {
      infoPtr = infoPtrIn; settingsPtr = settingsPtrIn; rndmPtr = rndmPtrIn;
      particleDataPtr = particleDataPtrIn;
      crossSectionDataPtr = crossSectionDataPtrIn; userHooksPtr = userHooksIn;
    }

  // @TODO: How to intialize decays?
  void init(Couplings* couplingsPtrIn, TimeShower* timesDecPtrIn, 
            DecayHandler* decayHandlerPtrIn, vector<int> handledParticles) {
    flavSel.init(*settingsPtr, particleDataPtr, rndmPtr, infoPtr);
    decays.init(infoPtr, *settingsPtr, particleDataPtr, rndmPtr, 
                couplingsPtrIn, timesDecPtrIn, &flavSel, decayHandlerPtrIn,
                handledParticles);
    resonanceDecays.init(infoPtr, particleDataPtr, rndmPtr);

    doSecondRescattering = settingsPtr->flag("Rescattering:doSecondRescattering");
    doDecays = settingsPtr->flag("Rescattering:doDecays");
    tau0Max = settingsPtr->parm("Rescattering:tau0Max");
    radiusMax = settingsPtr->parm("Rescattering:radiusMax");
  }

  void next(Event& event);

private: 

  struct PriorityVertex;

  Info* infoPtr;

  Settings* settingsPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  CrossSectionData* crossSectionDataPtr;

  UserHooks* userHooksPtr;

  StringFlav flavSel;

  ParticleDecays decays;

  ResonanceDecays resonanceDecays;
 
  void calcDecaysRescatters(Event& event, int iStart,
                            priority_queue<PriorityVertex>& queue);
  
  bool calculateInteraction(int idA, int idB, Event& event, Vec4& originOut);
  
  bool produceDecayProducts(int iDec, Event& event);

  void produceScatteringProducts(int iP1, int iP2, Vec4& origin, Event& event);

  bool doSecondRescattering, doDecays;
  double tau0Max, radiusMax;

};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
