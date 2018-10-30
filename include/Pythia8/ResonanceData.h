#ifndef Resonance_Data_H
#define Resonance_Data_H

#include "Pythia8/Interpolator.h"
#include "Pythia8/MassDependentWidth.h"
#include "Pythia8/Pythia.h"

namespace Pythia8 {


// @TODO temporary - to help thinking about the architecture
typedef string ResClass;
typedef string ResGenus;

class ResonanceData {
public:

  ResonanceData() {
    ifstream stream("ParticleWidths.xml");
    particleWidthPtr = new MassDependentWidth;
    particleWidthPtr->readXML(stream);
  }

  void initPtr(ParticleData* particleDataPtrIn) {
    particleDataPtr = particleDataPtrIn;
  }

  bool readXML(istream& inStream);


  const Interpolator& getDiffractiveSigmaDistribution(int idA, int idB) const;

  double getDiffractiveSigma(int idA, int idB, double eCM) const;

  double getBR(int idR, int idA, int idB) const;

  int getIsospin(int species) const;

  double getClebschGordan(int, int , int , int , int , int ) const ;
// @TODO Probably take Particles instead of just ids
  double getResonanceSigma(int idA, int idB, double eCM) const;

  vector<pair<pair<int, int>, double>> getOutputsWithFrequencies(int idA, int idB, double eCM) const {
    
    auto outs = getPossibleOutputs(idA, idB);
    vector<pair<pair<int, int>, double>> outsWithFreq(outs.size());

    for (int iR = 0; iR < (int)outs.size(); ++iR) {
      auto gens = genify(outs[iR].first, outs[iR].second);
      double weight = partialSigmaDistribution.at(gens)(eCM);
      outsWithFreq[iR] = pair<pair<int, int>, double>(outs[iR], weight);
    }

    return outsWithFreq;
  }

  vector<pair<int, int>> getPossibleOutputs(int idA, int idB) const;

  void print() const;

  void sanityCheck();

private:

  ParticleData* particleDataPtr;

  MassDependentWidth* particleWidthPtr;

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

  map<ResGenus, int> isospinType;

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