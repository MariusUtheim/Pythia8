
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergyResonance.h"


namespace Pythia8 {

//==========================================================================

// LowEnergyResonance class.
// This deals with cross sections and scattering through resonances.

//--------------------------------------------------------------------------

bool LowEnergyResonance::init(string path) {
  
  ifstream stream(path);
  if (!stream.is_open()) {
    infoPtr->errorMsg( "Warning in LowEnergyResonance::init: "
      "cannot open file");
    return false;
  }

  if (!particleWidths.readXML(stream)) {
    infoPtr->errorMsg( "Warning in LowEnergyResonance::init: "
      "cannot read resonance widths");
    return false;
  }

  for (auto id : particleWidths.getResonances()) {
    //@TODO: If a particle in particleWidths is not a hadron, it's an error
    if (!particleDataPtr->isHadron(id)
        || particleDataPtr->heaviestQuark(id) > 3)
      continue;
    
    // Signature takes the form BQS
    int charge = particleDataPtr->chargeType(id), strangeness = particleDataPtr->nStrangeQuarks(id);
    int signature = 100 * (particleDataPtr->isBaryon(id))
                  +  10 * ((charge >= 0) ? charge : (10 + charge))
                  +   1 * ((strangeness >= 0) ? strangeness : (10 + strangeness));

    auto iter = signatureToParticles.find(signature);
    if (iter != signatureToParticles.end())
      iter->second.push_back(id);
    else
      signatureToParticles.emplace(signature, vector<int>{id});
  }

  return true;
}

//--------------------------------------------------------------------------

int LowEnergyResonance::pickResonance(int idA, int idB, double eCM) {

  // Find possible resonances and their relative probabilities
  vector<int> candidates = getPossibleResonances(idA, idB);
  if (candidates.size() == 0) 
    return (cout << "Got no resonances" << endl), -1;

  vector<double> probs(candidates.size());
  for (size_t i = 0; i < candidates.size(); ++i)
    probs[i] = getPartialResonanceSigma(idA, idB, candidates[i], eCM);

  return candidates[rndmPtr->pick(probs)];
}

vector<int> LowEnergyResonance::pickDecayProducts(int idRes, double eCM) {

  auto brs = particleWidths.getWeightedProducts(idRes, eCM);
  if (brs.size() == 0) {
    // @TODO: This would be a bug
    cout << "Got no decay modes" << endl;
    return vector<int>();
  }

  vector<double> weights(brs.size());
  for (size_t i = 0; i < brs.size(); ++i)
    weights[i] = brs[i].first;

  return brs[rndmPtr->pick(weights)].second;
}

//--------------------------------------------------------------------------

double LowEnergyResonance::getPartialResonanceSigma(int idA, int idB, int idR, double eCM) const {
  
  // Get mass dependent width. If it is zero, the resonance cannot be formed
  double gammaR = particleWidths.width(idR, eCM);
  if (gammaR == 0)
    return 0.;

  // @TODO: Ordering matters, make sure the ordering of ids is canonical.
  double branchingRatio = particleWidths.branchingRatio(idR, vector<int>{idA, idB}, eCM);
  if (branchingRatio == 0)
    return 0.;

  double s = pow2(eCM), mR = particleDataPtr->m0(idR),
         mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);

  double pCMS2 = 1 / (4 * s) * (s - pow2(mA + mB)) * (s - pow2(mA - mB));

  // Calculate the resonance sigma
  double sigma = 
      M_PI / pCMS2
    * particleDataPtr->spinType(idR)
        / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB))
    * branchingRatio * pow2(gammaR) / (pow2(mR - eCM) + 0.25 * pow2(gammaR))
    * GEVINVSQ2MB;

  return sigma;
}

//--------------------------------------------------------------------------

double LowEnergyResonance::getResonanceSigma(int idA, int idB, double eCM) const {

  // For K_S and K_L, take average of K and Kbar
  if (idA == 310 || idA == 130)
    return 0.5 * (getResonanceSigma(311, idB, eCM) + getResonanceSigma(-311, idB, eCM));
  if (idB == 310 || idB == 130)
    return 0.5 * (getResonanceSigma(idA, 311, eCM) + getResonanceSigma(idA, -311, eCM));

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);

  // If the baryon is an antiparticle, flip both ids
  if (idA < 0) {
    idA = -idA;
    if (particleDataPtr->hasAnti(idB))
      idB = -idB;
  }
  
  // Abort if total energy is too low
  if (eCM < particleDataPtr->m0(idA) + particleDataPtr->m0(idB)) 
    return 0.;

  // Sum over all possible resonances
  vector<int> resonanceCandidates = getPossibleResonances(idA, idB);

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    // @TODO If we don't need partial sigma to be public, just put that code
    //       here instead
    sigmaRes += getPartialResonanceSigma(idA, idB, idR, eCM);
  }

  return sigmaRes;
}

//--------------------------------------------------------------------------

vector<int> LowEnergyResonance::getPossibleResonances(int idA, int idB) const {

  // Calculate total signature
  int baryonNumber = particleDataPtr->isBaryon(idA)
                   + particleDataPtr->isBaryon(idB); // @TODO: What about antiparticles?
  int charge = particleDataPtr->chargeType(idA) 
             + particleDataPtr->chargeType(idB);
  int strangeness = particleDataPtr->nStrangeQuarks(idA) 
                  + particleDataPtr->nStrangeQuarks(idB);

  int signature = 100 * (baryonNumber)
                +  10 * ((charge >= 0) ? charge : (10 + charge))
                +   1 * ((strangeness >= 0) ? strangeness : (10 + strangeness));

  // Get the resonances that conserve the signature
  auto candidates = signatureToParticles.find(signature);
  if (candidates == signatureToParticles.end())
    return vector<int>();
  else
    return candidates->second;
}

//==========================================================================

}