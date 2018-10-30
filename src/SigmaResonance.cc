
#include "Pythia8/SigmaResonance.h"

namespace Pythia8 { 

//SigmaResonance::SigmaResonance()
//  : sigmaNucleonTotal("sigmaTotal.dat"),
//    sigmaNN2NDelta("sigmaNN2NDelta.dat"),
//    sigmaNN2DeltaDelta("sigmaNN2DeltaDelta.dat"),
//    sigmaNN2NNstar("sigmaNN2NNstar.dat"),
//    sigmaNN2NDeltastar("sigmaNN2NDeltastar.dat"),
//    sigmaNN2DeltaNstar("sigmaNN2DeltaNstar.dat"),
//    sigmaNN2DeltaDeltastar("sigmaNN2DeltaDeltastar.dat")
//{ } 




double SigmaResonance::sigma(int idA, int idB, double eCM) const {
  auto& distribution = resonanceDataPtr->getDiffractiveSigmaDistribution(idA, idB);
  return distribution(eCM);
}

pair<int, int> SigmaResonance::pickProducts(int idA, int idB, double eCM) {
  
  auto candidates = resonanceDataPtr->getOutputsWithFrequencies(idA, idB, eCM);

  vector<double> weights(candidates.size());

  for (size_t i = 0; i < candidates.size(); ++i)
    weights[i] = candidates[i].second;

  int chosenIndex = rndmPtr->pick(weights);

  return candidates[chosenIndex].first;
}


}
