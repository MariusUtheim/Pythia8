
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


enum ResonanceClass {
  N = 1, D = 2, Nstar = 3, Dstar = 4,
  aN = -1, aD = -2, aNstar = -3, aDstar = -4
};


static ResonanceClass getClass(int id) {
  throw "Not implemented";
}

static vector<pair<int, int>> getScatteringResonances(pair<int, int> inputs) {
  /* 
    if resA > 0 && resB > 0 then
      resonances
    elif resA < 0 && resB < 0 then
      [ for res in resonances -> -res ]
    else
      [ for (resC, resD) in resonances do
        yield (-resC, resD)
        if resC <> resD then yield (resC, resD) ]
  */ 
  throw "Not implemented";
}

static double _sigma(ResonanceClass resA, ResonanceClass resB, ResonanceClass resC, ResonanceClass resD) {
  throw "Not implemented";
}

static vector<int> getParticlesFromClass(ResonanceClass cls) {
  throw "Not implemented";
}

static double sigmaFromTable(ResonanceClass resA, ResonanceClass resB, double eCM) {
  throw "Not implemented";
}
static double sigmaFromTable(ResonanceClass resA, ResonanceClass resB, ResonanceClass resC, ResonanceClass resD, double eCM) {
  throw "Not implemented";
}

pair<int, int> SigmaResonance::getResonanceClasses(int idA, int idB) const {
  
  ResonanceClass resA = getClass(idA), resB = getClass(idB);
  if (resA < resB) {
    swap(resA, resB);
  }
  return pair<int, int>(resA, resB);
}

double SigmaResonance::sigma(int idA, int idB, double eCM) const {
  pair<int, int> entry = getResonanceClasses(idA, idB);
  return sigmaTotals.at(entry)(eCM);
}

vector<int> SigmaResonance::pickProducts(int idAIn, int idBIn, double eCM) {
  pair<int, int> entry = getResonanceClasses(idAIn, idBIn);

  vector<pair<int, int>> resonances = getScatteringResonances(entry);
  vector<double> weights(resonances.size());

  for (int i = 0; i < resonances.size(); ++i) {
    weights[i] = brFactor.at(resonances[i])(eCM);
  }
}


}
