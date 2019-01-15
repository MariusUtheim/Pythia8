#ifndef Low_Energy_Resonance_H
#define Low_Energy_Resinance_H

#include "Pythia8/Event.h"
#include "Pythia8/LowEnergyData.h"
#include "Pythia8/LowEnergyProcess.h"

namespace Pythia8 {

class LowEnergyResonance: public LowEnergyProcess {
public:

  void initPtr(ParticleData* particleDataPtrIn, LowEnergyData* lowEnergyDataPtrIn)
  { particleDataPtr = particleDataPtrIn; lowEnergyDataPtr = lowEnergyDataPtrIn; }

  double sigmaTotal(int i1, int i2, const Event& event) const { return 0.; }

  bool collide(int i1, int i2, Event& event) { return false; }

  
  double getPartialResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;

  double getPartialElasticResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;

  double getResonanceSigma(int idA, int idB, double eCM) const;

  double getElasticResonanceSigma(int idA, int idB, double eCM) const;


private:

  ParticleData* particleDataPtr;

  LowEnergyData* lowEnergyDataPtr;

};

}

#endif