
#include "Pythia8/LowEnergyController.h"


namespace Pythia8 {

//--------------------------------------------------------------------------

double LowEnergyController::getTotalSigma(int idA, int idB, double eCM) const {

  if (eCM < particleDataPtr->m0(idA) + particleDataPtr->m0(idB))
    return 0.;

  if (particleDataPtr->isBaryon(idA) && particleDataPtr->isBaryon(idB)) {
    if (idA * idB < 0) {
      if (idA > 0)
        return getTotalSigmaBBbar(idA, idB, eCM);
      else
        return getTotalSigmaBBbar(idB, idA, eCM);
    }
    else
      return getTotalSigmaBB(idA, idB, eCM);
  }
  else {
    if (particleDataPtr->isMeson(idB))
      return getTotalSigmaXM(idA, idB, eCM);
    else
      return getTotalSigmaXM(idB, idA, eCM);
  }
  
  return 0.;
}


//--------------------------------------------------------------------------

// B+B section

double LowEnergyController::getTotalSigmaBB(int idA, int idB, double eCM) const {

  if (idA == 2212 && idB == 2212)
    return lowEnergyDiffractive.getDiffractiveSigma(idA, idB, eCM);

  cout << "LowEnergyController::getTotalSigmaBB not implemented" << endl;
  return 0.;
}

/*
const Interpolator& LowEnergyController::getDiffractiveSigmaDistribution(int idA, int idB) const {
  auto cls = classify(idA, idB);
  auto ptr = totalSigmaDistribution.find(cls);
  if (ptr == totalSigmaDistribution.end())
    return Interpolator::Zero;
  else
    return ptr->second;
}



vector<pair<int, int>> LowEnergyController::getDiffractiveOutputs(int idA, int idB) const {
  ResClass clsA = speciesToClass(idA), clsB = speciesToClass(idB);
  auto cls = pair<ResClass, ResClass>(clsA, clsB);

  int totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  auto& outputClasses = scatterChannels.at(cls);

  vector<pair<int, int>> possibleOutputs;

  for (auto outputClass : outputClasses) {

    auto parCs = classToSpecies(outputClass.first), 
         parDs = classToSpecies(outputClass.second);

    if (idA > 0 && idB > 0) {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge) {
            possibleOutputs.push_back(pair<int, int>(pC, pD));

          }
    }
    else if (idA < 0 && idB < 0) {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge)
            possibleOutputs.push_back(pair<int, int>(-pC, -pD));
    }
    else {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge) {
            possibleOutputs.push_back(pair<int, int>(-pC, pD));
            possibleOutputs.push_back(pair<int, int>(pC, -pD));
          }
    }

  }
  
  return possibleOutputs;
}

*/

//--------------------------------------------------------------------------

// B+Bbar section

double LowEnergyController::getTotalSigmaBBbar(int idA, int idB, double eCM) const {
  cout << "LowEnergyController::getTotalSigmaBBbar not implemented" << endl;
  return 0.;
}
/*
// Calculate sigma_annihilation. Assumes one index refers to a baryon, the other to an antibaryon
double LowEnergyController::getAnnihilationSigma(int idA, int idB, double eCM) const {
  if (idA > idB) 
    swap(idA, idB);

  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM <= mA + mB)
    return 0.;

  double xsA = getStrangenessFactor(idA), xsB = getStrangenessFactor(idB);

  static const double sigma0N = 120., s0 = 3.5214, A2 = 0.0025, B = 0.6;
  double s = eCM * eCM;
  double sigmaNN = sigma0N * s0/s * (A2 * s0 / (pow2(s - s0) + A2 * s0) + B);

  // @TODO separate into annihilation, diffractive and elastic

  return sigmaNN * (1. - 0.4 * xsA) * (1. - 0.4 * xsB);
}
*/
//--------------------------------------------------------------------------

static const vector<double> ppiminusInterpolData = 
  { 20.3919, 35.1056, 53.4896, 67.8302, 61.3948, 47.1339, 36.1912, 
    29.9881, 26.8335, 26.5702, 27.7574, 28.5599, 30.0137, 32.8038, 
    37.3226, 44.0601, 45.7762, 42.1849, 37.6042, 35.9693, 38.7317, 
    45.5506, 53.6462, 58.0236, 54.3861, 47.0179, 41.6299, 38.2877, 
    36.5429, 36.1408, 36.1086, 36.3178, 36.1149, 35.7958, 35.3629, 
    34.9317, 34.5016, 34.3947, 34.5755, 34.7563, 34.9371, 35.1311, 
    35.3531, 35.575, 35.797, 35.9509, 35.7146, 35.4783, 35.2419, 35.0056, 
    34.7693, 34.5329, 34.2954, 34.0514, 33.8073, 33.5632, 33.3192, 
    33.0751, 32.831, 32.587, 32.381, 32.278, 32.1749, 32.0718, 31.9688, 
    31.8657, 31.7626, 31.6595, 31.5565, 31.4534, 31.3503, 31.244, 
    31.1329, 31.0218, 30.9108, 30.7997, 30.6887, 30.5776, 30.4666, 
    30.3555, 30.2445, 30.1334, 30.0257, 29.9217, 29.8176, 29.7135, 
    29.6095, 29.5054, 29.4014, 29.2973, 29.1932, 29.0978, 29.0211, 
    28.9445, 28.8678, 28.7911, 28.7144, 28.6378, 28.5611, 28.4844, 
    28.4077, 28.3311, 28.2806, 28.2432, 28.2058, 28.1684, 28.131, 
    28.0936, 28.0562, 28.0188, 27.9814, 27.944, 27.9066, 27.8692, 
    27.8318, 27.7944, 27.757, 27.7196, 27.6822, 27.6448, 27.6074 };

static const Interpolator ppiminusInterpol(1.15, 3.89, ppiminusInterpolData);

static const vector<double> ppiplusInterpolData =
 { 33.1408, 30.3605, 29.3346, 29.8405, 30.06, 30.1206, 29.917, 29.3581,
   28.8239, 28.5466, 28.2694, 27.9941, 27.7218, 27.4494, 27.1523,
   26.8496, 26.547, 26.3317, 26.1572, 25.9828, 25.8343 };
static const Interpolator ppiplusInterpol(2.03, 3.50, ppiplusInterpolData);


double LowEnergyController::getTotalSigmaXM(int idX, int idM, double eCM) const {
  if (idX == 2212 && idM == -211) {
    return eCM <= 1.8 ? lowEnergyResonance.getResonanceSigma(idX, idM, eCM)
         : eCM <= 3.5 ? ppiminusInterpol(eCM)
         :              13.7 * pow(eCM, 0.158) + 35.9 * pow(eCM, -0.90);
  }
  else if (idX == 2212 && idM == 211) {
    return eCM <= 2.03 ? lowEnergyResonance.getResonanceSigma(idX, idM, eCM)
         : eCM <= 3.50 ? ppiplusInterpol(eCM)
         :              13.7 * pow(eCM, 0.158) + 35.9 * pow(eCM, -0.90);
  }
  else
    return 0.;
}

const LowEnergyProcess& LowEnergyController::pickProcess(int idA, int idB, double eCM) {

  if (particleDataPtr->isMeson(idB)) {
    double sigmaTotal = getTotalSigmaXM(idA, idB, eCM);
    double pRes = lowEnergyResonance.getResonanceSigma(idA, idB, eCM) / sigmaTotal;
    double pEl = lowEnergyResonance.getElasticResonanceSigma(idA, idB, eCM) / sigmaTotal;
    double pStr = 1 - pRes - pEl;

    switch (rndmPtr->pick({ pRes, pEl, pStr })) {
      case 0: return lowEnergyResonance;
      case 1: return lowEnergyResonance;// @TODO return lowEnergyElastic;
      case 2: return lowEnergyStrings;
      default: throw "Error in Pythia internal logic (LowEnergyController)"; // @TODO
    }
  }
  else if (idA == 2212 && idB == 2212) {
    return lowEnergyDiffractive;
  }
  else {
    cout << "NOT IMPLEMENTED:" << __LINE__ << endl;
  }
}


//--------------------------------------------------------------------------


}
