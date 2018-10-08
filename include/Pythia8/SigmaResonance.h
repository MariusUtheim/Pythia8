#ifndef Pythia8_SigmaResonance_H
#define Pythia8_SigmaResonance_H

#include "Pythia8/Interpolator.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

class SigmaResonance {

public:

  SigmaResonance();

  double sigma(int idA, int idB, double eCM) const;

  vector<int> pickProducts(int idAIn, int idBIn, double eCM);


  const Interpolator& sigmaDistribution(int idA, int idB) const; 
  const Interpolator& sigmaDistribution(int idA, int idB, int idC, int idD) const;

private:

  ParticleData* particleDataPtr;

  map<pair<int, int>, Interpolator&> sigmaTotals;

  map<pair<int, int>, Interpolator&> brFactor;

  pair<int, int> getResonanceClasses(int idA, int idB) const;

};

}

#endif