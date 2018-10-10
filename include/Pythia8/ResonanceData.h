#ifndef Resonance_Data_H
#define Resonance_Data_H

#include "Pythia8/Interpolator.h"
#include "Pythia8/Pythia.h"

namespace Pythia8 {


// @TODO temporary - to help thinking about the architecture
typedef string ResClass;
typedef string ResGenus;

class ResonanceData {
public:

  ResonanceData() {}

  void initPtr(ParticleData* particleDataPtrIn) {
    particleDataPtr = particleDataPtrIn;
  }

  bool readXML(istream& inStream);


  double getTotalSigma(int idA, int idB, double eCM) const {
    auto cls = classify(idA, idB);
    auto ptr = totalSigmaDistribution.find(cls);
    if (ptr == totalSigmaDistribution.end())
      return 0.;
    else
      return ptr->second(eCM);
  }

  vector<pair<pair<int, int>, double>> getOutputsWithFrequencies(int idA, int idB, double eCM) const {
    auto outs = getPossibleOutputs(idA, idB);

    vector<pair<pair<int, int>, double>> outsWithFreq(outs.size());

    for (int i = 0; i < (int)outs.size(); ++i) {
      auto gens = genify(outs[i].first, outs[i].second);
      outsWithFreq[i] = pair<pair<int, int>, double>(outs[i], partialSigmaDistribution.at(gens)(eCM));
    }

    return outsWithFreq;
  }

  const Interpolator& getTotalSigmaDistribution(int idA, int idB) const {
    auto cls = classify(idA, idB);
    auto ptr = totalSigmaDistribution.find(cls);
    if (ptr == totalSigmaDistribution.end())
      return Interpolator::Zero;
    else
      return ptr->second;
  }

  vector<pair<int, int>> getPossibleOutputs(int idA, int idB) const;

  void print() const;

  void sanityCheck();

private:

  ParticleData* particleDataPtr;

  ResClass classify(int id) const {
    return genusToClass.at(speciesToGenus.at(abs(id)));
  }

  ResGenus genify(int id) const {
    return speciesToGenus.at(abs(id));
  }

  pair<ResClass, ResClass> classify(int idA, int idB) const {
    ResClass first = classify(idA);
    ResClass second = classify(idB);
    return pair<ResClass, ResClass>(first, second);
  }

  pair<ResGenus, ResGenus> genify(int idA, int idB) const {
    ResClass first = genify(idA);
    ResClass second = genify(idB);
    return pair<ResGenus, ResGenus>(first, second);
  }

  map<pair<ResClass, ResClass>, vector<pair<ResClass, ResClass>>> scatterChannels;

  map<pair<ResClass, ResClass>, Interpolator> totalSigmaDistribution;

  map<pair<ResGenus, ResGenus>, Interpolator> partialSigmaDistribution;

  map<ResGenus, ResClass> genusToClass;
  map<int, ResGenus> speciesToGenus;
  map<ResClass, vector<ResGenus>> classToGenera;
  map<ResGenus, vector<int>> genusToParticles;

  vector<int> classToSpecies(ResClass cls) const {
    auto genera = classToGenera.at(cls);
    vector<int> species;
    for (auto gen : genera)
      for (auto spc : genusToParticles.at(gen))
        species.push_back(spc);
    return species;
  }
};

}

#endif