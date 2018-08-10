
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
  Hist eCM;

  static RescatteringLogger& _()
  {
    static RescatteringLogger _instance;
    return _instance;
  }

  RescatteringLogger(RescatteringLogger const &) = delete;
  void operator =(RescatteringLogger const &) = delete;

  void print()
  {
    int nTot = nHit + nObviousAway + nMissImpact + nMissAway + nSpacelike;
    double totFac = sqrt(nHit) + sqrt(nObviousAway) + sqrt(nMissImpact) + sqrt(nMissAway) + sqrt(nSpacelike);
    cout << "\n ==== Rescattering report === " << endl
         << left 
         << " #nHit   = " << setw(10) << nHit        << 100 * sqrt(nHit)        / totFac << endl
         << " #Away!  = " << setw(10) << nObviousAway<< 100 * sqrt(nObviousAway)/ totFac << endl
         << " #TooFar = " << setw(10) << nMissImpact << 100 * sqrt(nMissImpact) / totFac << endl
         << " #Away   = " << setw(10) << nMissAway   << 100 * sqrt(nMissAway  ) / totFac << endl
         << " #S.like = " << setw(10) << nSpacelike  << 100 * sqrt(nSpacelike ) / totFac << endl
         << " #TOTAL  = " << setw(1) << (int)nTot << endl
         << impactParameter
         << r << phi
         << eCM;


    HistPlot hpl("myplot");
    hpl.plotFrame("outplot", impactParameter / nHit, "Impact parameter", "$b$", "$p$");
    hpl.plotFrame("", r / nHit, "Collision radius", "$r$", "$p$");
    hpl.plotFrame("", eCM / nHit, "Center of mass energy", "$E_{CM}$", "$p$");
  }

private:
  RescatteringLogger() 
    : nMissImpact(0), nMissAway(0), nSpacelike(0),
      impactParameter("Impact parameter", 50, 0, 4, false),
      r("r", 50, 0, 100, false),
      phi("phi", 100, -M_PI, M_PI),
      eCM("eCM", 30, 0.05, 5, true)
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
  // @TODO: Cache settings
	void init(Couplings* couplingsPtrIn, TimeShower* timesDecPtrIn, 
					  DecayHandler* decayHandlerPtrIn, vector<int> handledParticles) {
		flavSel.init(*settingsPtr, particleDataPtr, rndmPtr, infoPtr);
		decays.init(infoPtr, *settingsPtr, particleDataPtr, rndmPtr, 
								couplingsPtrIn, timesDecPtrIn, &flavSel, decayHandlerPtrIn,
								handledParticles);
		resonanceDecays.init(infoPtr, particleDataPtr, rndmPtr);
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

};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
