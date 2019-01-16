#ifndef Low_Energy_Controller_H
#define Low_Energy_Controller_H

#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergyDiffractive.h"
#include "Pythia8/LowEnergyResonance.h"
#include "Pythia8/MassDependentWidth.h"

namespace Pythia8 {


// @TODO temporary - to help thinking about the architecture
typedef string ResClass;
typedef string ResGenus;

class LowEnergyController {
public:

  void initPtr(Rndm* rndmPtrIn, ParticleData* particleDataPtrIn, LowEnergyData* lowEnergyDataPtrIn)
  {
    rndmPtr = rndmPtrIn;
    particleDataPtr = particleDataPtrIn;
    lowEnergyResonance.initPtr(rndmPtrIn, particleDataPtrIn, lowEnergyDataPtrIn);
    lowEnergyDiffractive.initPtr(rndmPtrIn, particleDataPtrIn, lowEnergyDataPtrIn);
  }


  double getTotalSigma(int idA, int idB, double eCM) const;

/*
  const Interpolator& getDiffractiveSigmaDistribution(int idA, int idB) const;

  double getAnnihilationSigma(int idA, int idB, double eCM) const;

  double getElasticSigma(int idA, int idB, double eCM) const;

*/


  const LowEnergyProcess& pickProcess(int idA, int idB, double eCM);

  void showPickProbabilities(int idA, int idB, double eCM) const;


private:

  double getTotalSigmaBB(int idA, int idB, double eCM) const;
  double getTotalSigmaBBbar(int idB, int idBbar, double eCM) const;
  double getTotalSigmaXM(int idX, int idM, double eCM) const;

  Rndm* rndmPtr; 
  
  ParticleData* particleDataPtr;

  LowEnergyResonance lowEnergyResonance;
  LowEnergyDiffractive lowEnergyDiffractive;
};

}

#endif