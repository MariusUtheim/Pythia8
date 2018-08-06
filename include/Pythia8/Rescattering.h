
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
	bool calculateInteraction(Particle& p1In, Particle& p2In, Vec4& originOut);
	
	bool produceDecayProducts(int iDec, Event& event);

	void produceScatteringProducts(int iP1, int iP2, Vec4& origin, Event& event);

};

} // end namespace Pythia8

#endif // Pythia8_Rescattering_H
