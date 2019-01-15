#ifndef Low_Energy_Controller_H
#define Low_Energy_Controller_H

#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergyResonance.h"
#include "Pythia8/MassDependentWidth.h"

namespace Pythia8 {


// @TODO temporary - to help thinking about the architecture
typedef string ResClass;
typedef string ResGenus;

class LowEnergyController {
public:

  void initPtr(ParticleData* particleDataPtrIn, LowEnergyData* lowEnergyDataPtrIn)
  {
    particleDataPtr = particleDataPtrIn;
    lowEnergyResonance.initPtr(particleDataPtrIn, lowEnergyDataPtrIn);
  }


  double getTotalSigma(int idA, int idB, double eCM) const;

/*
  const Interpolator& getDiffractiveSigmaDistribution(int idA, int idB) const;

  double getDiffractiveSigma(int idA, int idB, double eCM) const;

  vector<pair<pair<int, int>, double>> getOutputsWithFrequencies(int idA, int idB, double eCM) const;

  vector<pair<int, int>> getDiffractiveOutputs(int idA, int idB) const;


  double getAnnihilationSigma(int idA, int idB, double eCM) const;

  double getElasticSigma(int idA, int idB, double eCM) const;

*/


  vector<Particle> pickProducts(int idA, int idB, double eCM) const;

  void showPickProbabilities(int idA, int idB, double eCM) const;


private:

  double getTotalSigmaBB(int idA, int idB, double eCM) const;
  double getTotalSigmaBBbar(int idB, int idBbar, double eCM) const;
  double getTotalSigmaXM(int idX, int idM, double eCM) const;

  ParticleData* particleDataPtr;

  LowEnergyResonance lowEnergyResonance;
};

}

#endif