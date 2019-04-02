
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
    if (!particleDataPtr->isHadron(id) || hasHeavyQuark(id))
      continue;

    // Signature takes the form BQS
    int charge = particleDataPtr->charge(id), strangeness = getStrangeness(id);
    int signature = (particleDataPtr->isBaryon(id) * 100)
                  + ((2 * abs(charge) - (charge < 0)) * 10)
                  + ((2 * abs(strangeness) - (strangeness < 0)));

    auto ptr = signatureToParticles.lower_bound(signature);
    if (ptr != signatureToParticles.end() && ptr->first == signature)
      ptr->second.push_back(id);
    else
      signatureToParticles.insert(ptr, make_pair(signature, vector<int>{id}));
  }

  return true;
}

//--------------------------------------------------------------------------

bool LowEnergyResonance::collide(int i1, int i2, Event& event, Vec4 origin) {
  Particle& hA = event[i1];
  Particle& hB = event[i2];

  double eCM = (hA.p() + hB.p()).mCalc();

  // Find possible resonances and their relative probabilities
  vector<int> candidates = getResonanceCandidates(hA.id(), hB.id());
  if (candidates.size() == 0) 
    return (cout << "Got no resonances" << endl), false;

  vector<double> probabilities(candidates.size());

  for (size_t i = 0; i < candidates.size(); ++i)
    (probabilities[i] = getPartialResonanceSigma(hA.id(), hB.id(), candidates[i], eCM)), (cout << probabilities[i] << " ");

  // Create the resonance
  int iNew = event.append(candidates[rndmPtr->pick(probabilities)],
                          115, i1, i2, 0, 0, 0, 0, hA.p() + hB.p(), eCM);
  event[iNew].vProd(origin);

  for (int i : { i1, i2 }) {
    event[i].daughters(iNew, iNew);
    event[i].statusNeg();
  }

  // Decay the resonance
  // @TODO: Do this properly
  auto brs = particleWidths.getWeightedProducts(event[iNew].id(), eCM);

  if (brs.size() == 0) {
    cout << "Got no decay modes" << endl;
    return false;
  }
  for (auto prpr : brs) {
    for (auto prod : prpr.second)
      cout << prod << " ";
    cout << ": " << prpr.first << endl;
  }

  double threshold = rndmPtr->flat();
  double cumulativeSum = 0.;
  for (auto br : brs) {
    cumulativeSum += br.first;
    if (cumulativeSum >= threshold) {
      for (int id : br.second)
        event.append(id, 116, iNew, iNew, 0, 0, 0, 0, Vec4());
      break;
    }
  }

  return true;
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

  // Calculate the resonance sigma
  double sigma = 
      4 * M_PI * s / ((s - pow2(mA + mB)) * (s - pow2(mA - mB))) 
    * particleDataPtr->spinType(idR) / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB))
    * branchingRatio * pow2(gammaR) / (pow2(mR - eCM) + 0.25 * pow2(gammaR))
    * INVGEVSQ2MB;

  // If the two particles are the same except for I3, then multiply by 2
  // (e.g. for pi+pi- --> rho)
  // @TODO: Check that this test is correct for all cases. What about charm?
  // @TODO: Is this necessary? Shouldn't it be included in gamma already?
  int quarksA = (idA / 10) % 1000, quarksB = (idB / 10) % 1000;
  if (quarksA != quarksB && (idA - 10 * quarksA) == (idB - 10 * quarksB)
      && getStrangeness(idA) == getStrangeness(idB))
  {
    sigma *= 2;
  }

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
  vector<int> resonanceCandidates = getResonanceCandidates(idA, idB);

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    // @TODO If we don't need partial sigma to be public, just put that code
    //       here instead
    sigmaRes += getPartialResonanceSigma(idA, idB, idR, eCM);
  }

  return sigmaRes;
}

//--------------------------------------------------------------------------

vector<int> LowEnergyResonance::getResonanceCandidates(int idA, int idB) const {

  // Calculate total signature
  int baryonNumber = particleDataPtr->isBaryon(idA)
                   + particleDataPtr->isBaryon(idB); // @TODO: What about antiparticles?
  int charge = particleDataPtr->charge(idA) + particleDataPtr->charge(idB);
  int strangeness = getStrangeness(idA) + getStrangeness(idB);

  int signature = (baryonNumber * 100) 
                + ((2 * abs(charge) - (charge < 0)) * 10)
                + ((2 * abs(strangeness) - (strangeness < 0)));

  // Get the resonances that conserve the signature
  auto candidates = signatureToParticles.find(signature);
  if (candidates == signatureToParticles.end())
    return vector<int>();
  else
    return candidates->second;
}

//--------------------------------------------------------------------------

int LowEnergyResonance::getStrangeness(int id) const {
  // In baryon id xxxabcx, count number of occurrences of '3' in abc
  id = abs(id);
  int count = 0;
  for (int n = (id / 10) % 1000; n > 0; n /= 10) {
    int j = n % 10;
  
    if (j == 3)
      count++;
  }

  return count;
}

bool LowEnergyResonance::hasHeavyQuark(int id) const {
  for (int n = (abs(id) / 10) % 1000; n > 0; n /= 10)
    if (n % 10 > 3)
      return true;
  return false;
}

//==========================================================================

}