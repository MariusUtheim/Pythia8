
#include "Pythia8/LowEnergyDiffractive.h"

bool LowEnergyDiffractive::collide(int i1, int i2, Event& event) const {
  Particle& hA = event[i1];
  Particle& hB = event[i2];

  if (hA.id() == 2212 && hB.id() == 2212) {
    auto outputs = lowEnergyDataPtr->getOutputsWithFrequencies(hA.id(), hB.id(), (hA.p() + hB.p()).mCalc());
    double threshold = rndmPtr->flat();
    double sum = 0.;
    for (auto output : outputs) {
      cout << output.first.first << " + " << output.first.second << " : " << output.second << endl;
      sum += output.second;
      if (sum >= threshold) {
        for (int iNew : {output.first.first, output.first.second}) 
          event.append(iNew, 117, i1, i2, 0, 0, 0, 0, Vec4());
        break;
      }
    }
  }
  else {
    cout << "NOT IMPLEMENTED:" << __FILE__ << ":" << __LINE__ << endl;
    return false;
  }
}


double LowEnergyDiffractive::getDiffractiveSigma(int idA, int idB, double eCM) const {
  return lowEnergyDataPtr->getDiffractiveSigma(idA, idB, eCM);
}
