
#ifndef Pythia8_Rescattering_H
#define Pythia8_Rescattering_H

#include "Pythia8/Basics.h"
#include "Pythia8/CrossSectionData.h"
#include "Pythia8/Event.h"
#include "Pythia8/ParticleDecays.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/ResonanceDecays.h"

namespace Pythia8 {

class RescatteringLogger
{
public:

  int nHit,  nObviousAway, nMissImpact, nMissAway, nSpacelike;
  Hist impactParameter;
  Hist r, phi;
  Hist eCM, ytau;
  int nScatter;

  static RescatteringLogger& instance()
  {
    static RescatteringLogger _instance;
    return _instance;
  }

  RescatteringLogger(RescatteringLogger const &) = delete;
  void operator =(RescatteringLogger const &) = delete;

  void print()
  {
    //int nTot = nHit + nObviousAway + nMissImpact + nMissAway + nSpacelike;
    //double totFac = sqrt(nHit) + sqrt(nObviousAway) + sqrt(nMissImpact) + sqrt(nMissAway) + sqrt(nSpacelike);
    cout 
         << phi
         << r
         << eCM
         << ytau;


    HistPlot hpl("myplot");
    hpl.plotFrame("outplot", phi / nScatter, "Azimuthal angle", "$\\phi$", "$n$");
    hpl.plotFrame("", r / nScatter, "Collision radius", "$r$", "$p$");
    hpl.plotFrame("", eCM / nScatter, "Center of mass energy", "$E_{CM}$", "$p$");
    hpl.plotFrame("", ytau / nScatter, "ytau", "$y_\\tau$", "$p$");
  }

private:
  RescatteringLogger() 
    : nMissImpact(0), nMissAway(0), nSpacelike(0),
      impactParameter("Impact parameter", 50, 0, 4, false),
      r("r", 50, 0, 100, false),
      phi("phi", 30, -M_PI, M_PI),
      eCM("eCM", 30, 0.1, 10, true),
      ytau("ytau", 50, -10, 10, false),
      nScatter(0)
  { }

};

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
    tau0Max = settingsPtr->parm("Rescattering:tau0Max");
    radiusMax = settingsPtr->parm("Rescattering:radiusMax");
  }

  void next(Event& event);

private: 

  Info* infoPtr;

  Settings* settingsPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  CrossSectionData* crossSectionDataPtr;

  UserHooks* userHooksPtr;

  StringFlav flavSel;

  ParticleDecays decays;

  ResonanceDecays resonanceDecays;
 

  bool calculateDecay(Particle& pIn, Vec4& originOut);
  bool calculateInteraction(int idA, int idB, Event& event, Vec4& originOut);
  
  bool produceDecayProducts(int iDec, Event& event);

  void produceScatteringProducts(int iP1, int iP2, Vec4& origin, Event& event);

  bool canScatter(Particle& particle);

  bool doSecondRescattering;
  double tau0Max, radiusMax;

};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
