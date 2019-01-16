
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergyResonance.h"


namespace Pythia8 {



/*
void LowEnergyController::showPickProbabilities(int idX, int idM, double eCM) const {
  if (idX == 2212 && idM == -211) {
    double total = getTotalSigmaXM(idX, idM, eCM);
    
    double resonant = getResonanceSigma(idX, idM, eCM);
    
    double mp = particleDataPtr->m0(2212), mpi = particleDataPtr->m0(211);
    double pLab = sqrt((pow2(eCM) - pow2(mp + mpi)) * (pow2(eCM) - (mp - mpi)) / (2 * mp));
    double elastic = 1.76 + 11.2 * pow(pLab, -0.64) + 0.043 * pow2(log(pLab))
                   - getElasticResonanceSigma(idX, idM, eCM);

    double string = total - resonant - elastic;

    cout << setw(8) << "Total: " << total << endl
         << setw(8) << "Res: " << resonant << endl
         << setw(8) << "El: " << elastic << endl
         << setw(8) << "Str: " << string << endl
         << endl;
  
  }
}
*/

bool LowEnergyResonance::collide(int i1, int i2, Event& event) const {
  Particle& hA = event[i1];
  Particle& hB = event[i2];
  
  double eCM = (hA.p() + hB.p()).mCalc();

  vector<int> candidates = lowEnergyDataPtr->getResonanceCandidates(hA.id(), hB.id());
  vector<double> probabilities(candidates.size());

  for (size_t i = 0; i < candidates.size(); ++i)
    probabilities[i] = getPartialResonanceSigma(hA.id(), hB.id(), candidates[i], false, eCM);
  
  int iNew = event.append(candidates[rndmPtr->pick(probabilities)],
                          115, i1, i2, 0, 0, 0, 0, hA.p() + hB.p(), eCM);

  for (int i : { i1, i2 }) {
    event[i].daughters(iNew, iNew);
    event[i].statusNeg();
  }

  // @TODO Decay the resonance
  auto brs = lowEnergyDataPtr->massDependentBRs(event[iNew].id(), eCM);
  double threshold = rndmPtr->flat();
  double cumulativeSum = 0.;
  for (auto br : brs) {
    cumulativeSum += br.second;
    if (cumulativeSum >= threshold) {
      for (int id : br.first)
        event.append(id, 116, iNew, iNew, 0, 0, 0, 0, Vec4());
      break;
    }
  }

  return true;
}


double LowEnergyResonance::getPartialResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const {
  double br = lowEnergyDataPtr->getBR(idR, idA, idB, eCM);
  if (br == 0.)
    return 0.;

  double gammaRes2 = pow2(lowEnergyDataPtr->massDependentWidth(idR, eCM));

  int iA = lowEnergyDataPtr->getIsospin(idA);
  int i3A = lowEnergyDataPtr->getIso3(idA);

  int iB = lowEnergyDataPtr->getIsospin(idB);
  int i3B = lowEnergyDataPtr->getIso3(idB);
  
  int iR = lowEnergyDataPtr->getIsospin(idR);
  int i3R = lowEnergyDataPtr->getIso3(idR);

  double cg2 = lowEnergyDataPtr->getClebschGordan2(iA, i3A, iB, i3B, iR, i3R);

  if (isnan(cg2)) {
    cout << " for " << particleDataPtr->name(idA) << " + " << particleDataPtr->name(idB) << " --> " << particleDataPtr->name(idR) << endl
         << "     < " << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " > " << endl;
  }

  int nJRes = particleDataPtr->spinType(idR);
  double m0 = particleDataPtr->m0(idR);

  double contribution = cg2 * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);

  if (gensEqual && idA != idB)
    contribution *= 2;

  return contribution;
}

double LowEnergyResonance::getPartialElasticResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const {
  double br = lowEnergyDataPtr->getBR(idR, idA, idB, eCM);
  if (br == 0.)
    return 0.;

  double gammaRes2 = pow2(lowEnergyDataPtr->massDependentWidth(idR, eCM));

  int iA = lowEnergyDataPtr->getIsospin(idA);
  int i3A = lowEnergyDataPtr->getIso3(idA);

  int iB = lowEnergyDataPtr->getIsospin(idB);
  int i3B = lowEnergyDataPtr->getIso3(idB);
  
  int iR = lowEnergyDataPtr->getIsospin(idR);
  int i3R = lowEnergyDataPtr->getIso3(idR);

  double cg2 = lowEnergyDataPtr->getClebschGordan2(iA, i3A, iB, i3B, iR, i3R);

  if (isnan(cg2)) {
    cout << " for " << particleDataPtr->name(idA) << " + " << particleDataPtr->name(idB) << " --> " << particleDataPtr->name(idR) << endl
         << "     < " << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " > " << endl;
  }

  int nJRes = particleDataPtr->spinType(idR);
  double m0 = particleDataPtr->m0(idR);

  double contribution = cg2 * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);

  if (gensEqual && idA != idB)
    contribution *= 2;

  contribution *= br * cg2; // Elastic correction

  return contribution;
}

// @TODO Probably take Particles instead of just ids
double LowEnergyResonance::getResonanceSigma(int idA, int idB, double eCM) const {

  // For K_S and K_L, take average of K and Kbar
  if (idA == 310 || idA == 130)
    return 0.5 * (getResonanceSigma(311, idB, eCM) + getResonanceSigma(-311, idB, eCM));
  if (idB == 310 || idB == 130)
    return 0.5 * (getResonanceSigma(idA, 311, eCM) + getResonanceSigma(idA, -311, eCM));

  // @TODO Deal with what happens when there is an antibaryon
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);


  vector<int> resonanceCandidates = lowEnergyDataPtr->getResonanceCandidates(idA, idB);

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    sigmaRes += getPartialResonanceSigma(idA, idB, idR, lowEnergyDataPtr->gensEqual(idA, idB), eCM);
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  // @TODO define constant 0.38937966 = GeV^-2 to mb
  double preFactor = 0.38937966 * M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  sigmaRes *= preFactor;

  return sigmaRes;
}


double LowEnergyResonance::getElasticResonanceSigma(int idA, int idB, double eCM) const {

  // For K_S and K_L, take average of K and Kbar
  if (idA == 310 || idA == 130)
    return 0.5 * (getResonanceSigma(311, idB, eCM) + getResonanceSigma(-311, idB, eCM));
  if (idB == 310 || idB == 130)
    return 0.5 * (getResonanceSigma(idA, 311, eCM) + getResonanceSigma(idA, -311, eCM));

  // @TODO Deal with what happens when there is an antibaryon
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);


  vector<int> resonanceCandidates = lowEnergyDataPtr->getResonanceCandidates(idA, idB);

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    sigmaRes += getPartialElasticResonanceSigma(idA, idB, idR, lowEnergyDataPtr->gensEqual(idA, idB), eCM);
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  // @TODO define constant 0.38937966 = GeV^-2 to mb
  double preFactor = 0.38937966 * M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  sigmaRes *= preFactor;

  return sigmaRes;
}


}