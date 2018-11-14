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

  double getStrangeness(int id) const;

  double getTotalSigma(int idA, int idB, double eCM) const;


  const Interpolator& getDiffractiveSigmaDistribution(int idA, int idB) const;

  double getDiffractiveSigma(int idA, int idB, double eCM) const;


  double getBR(int idR, int idA, int idB, double eCM) const;

  int getIsospin(int species) const;

  int getIso3(int species) const;

  double getClebschGordan2(int, int , int , int , int , int ) const ;

  double getResonanceSigma(int idA, int idB, double eCM) const;

  double getStrangenessFactor(int id) const;

  double getAnnihilationSigma(int idA, int idB, double eCM) const;

  double getElasticSigma(int idA, int idB, double eCM) const;


  vector<pair<pair<int, int>, double>> getOutputsWithFrequencies(int idA, int idB, double eCM) const;

  vector<pair<int, int>> getDiffractiveOutputs(int idA, int idB) const;


  void print() const;

  bool sanityCheck();

private:

  ParticleData* particleDataPtr;

  MassDependentWidth* particleWidthPtr;

  ResClass classify(int id) const {
    return genusToClass.at(speciesToGenus.at(abs(id)));
  }

  ResGenus genify(int id) const {
    auto ptr = speciesToGenus.find(abs(id));
    if (ptr == speciesToGenus.end())
      return (cout << "Genus not found for particle species " << id << endl), "";
    else
      return ptr->second;
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