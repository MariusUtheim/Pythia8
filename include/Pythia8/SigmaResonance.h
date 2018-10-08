#ifndef Pythia8_SigmaResonance_H
#define Pythia8_SigmaResonance_H

#include "Pythia8/Interpolater.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

class SigmaResonance {

public:

  SigmaResonance();

  double sigma(int idA, int idB, double eCM) const;

  vector<int> pickProducts(int idAIn, int idBIn, double eCM);


  const Interpolater& sigmaDistribution(int idA, int idB) const; 
  const Interpolater& sigmaDistribution(int idA, int idB, int idC, int idD) const;

private:

  ParticleData* particleDataPtr;

  map<pair<int, int>, Interpolater&> sigmaTotals;

  map<pair<int, int>, Interpolater&> brFactor;

  pair<int, int> getResonanceClasses(int idA, int idB) const;

};

}

#endif