#ifndef Resonance_Data_H
#define Resonance_Data_H

#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/MassDependentWidth.h"

namespace Pythia8 {


// @TODO temporary - to help thinking about the architecture
typedef string ResClass;
typedef string ResGenus;

class ResonanceData {
public:

  bool init(string path) {
    ifstream stream(path);
    if (!stream.is_open()) return false;
    return readXML(stream);
  }

  void initPtr(ParticleData* particleDataPtrIn, MassDependentWidth* particleWidthPtrIn) {
    particleDataPtr = particleDataPtrIn; particleWidthPtr = particleWidthPtrIn;
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

  vector<int> getResonanceCandidates(int idA, int idB) const;
  double getPartialResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;
  double getPartialElasticResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const;
  double getResonanceSigma(int idA, int idB, double eCM) const;
  double getElasticResonanceSigma(int idA, int idB, double eCM) const;

  double getStrangenessFactor(int id) const;

  double getAnnihilationSigma(int idA, int idB, double eCM) const;

  double getElasticSigma(int idA, int idB, double eCM) const;


  vector<pair<pair<int, int>, double>> getOutputsWithFrequencies(int idA, int idB, double eCM) const;

  vector<pair<int, int>> getDiffractiveOutputs(int idA, int idB) const;

  vector<Particle> pickProducts(int idA, int idB, double eCM) const;

  void showPickProbabilities(int idA, int idB, double eCM) const;

  void print() const;

  bool sanityCheck();

private:

  double getTotalSigmaBB(int idA, int idB, double eCM) const;
  double getTotalSigmaBBbar(int idB, int idBbar, double eCM) const;
  double getTotalSigmaXM(int idX, int idM, double eCM) const;

  ParticleData* particleDataPtr;

  MassDependentWidth* particleWidthPtr;


  pair<ResClass, ResClass> classify(int idA, int idB) const {
    ResClass first = speciesToClass(idA);
    ResClass second = speciesToClass(idB);
    return pair<ResClass, ResClass>(first, second);
  }

  pair<ResGenus, ResGenus> genify(int idA, int idB) const {
    ResClass first = speciesToGenus(idA);
    ResClass second = speciesToGenus(idB);
    return pair<ResGenus, ResGenus>(first, second);
  }

  map<pair<ResClass, ResClass>, vector<pair<ResClass, ResClass>>> scatterChannels;

  map<pair<ResClass, ResClass>, Interpolator> totalSigmaDistribution;

  map<pair<ResGenus, ResGenus>, Interpolator> partialSigmaDistribution;

  map<ResGenus, ResClass> _genusToClass;
  ResClass genusToClass(ResGenus gen) const {
    auto ptr = _genusToClass.find(gen);
    if (ptr == _genusToClass.end())
      return (cout << "Class not found for particle genus " << gen << endl), "";
    else
      return ptr->second;
  }

  map<int, ResGenus> _speciesToGenus;
  ResGenus speciesToGenus(int spc) const {
    auto ptr = _speciesToGenus.find(abs(spc));
    if (ptr == _speciesToGenus.end())
      return (cout << "Genus not found for particle species " << spc << endl), "";
    else
      return ptr->second;
  }

  ResClass speciesToClass(int spc) const {
    auto gen = speciesToGenus(spc);
    if (gen == "") return "";
    return genusToClass(gen);
  }


  map<ResGenus, vector<int>> _genusToSpecies;
  vector<int> genusToSpecies(ResGenus gen) const {
    auto ptr = _genusToSpecies.find(gen);
    if (ptr == _genusToSpecies.end())
      return (cout << "Species not found for particle genus " << gen << endl), vector<int>();
    else
      return ptr->second;
  }

  map<ResClass, vector<ResGenus>> _classToGenera;
  vector<ResGenus> classToGenera(ResClass cls) const {
    auto ptr = _classToGenera.find(cls);
    if (ptr == _classToGenera.end())
      return (cout << "Genera not found for particle class " << cls << endl), vector<ResGenus>();
    else
      return ptr->second;
  }

  vector<int> classToSpecies(ResClass cls) const {
    auto genera = _classToGenera.find(cls);
    if (genera == _classToGenera.end())
      return (cout << "Genera not found for particle class " << cls << endl), vector<int>();
    
    vector<int> species;
    for (auto gen : genera->second) {
      auto spcs = _genusToSpecies.find(gen);
      if (spcs == _genusToSpecies.end())
        return (cout << "Species not found for particle genus " << gen << endl), vector<int>();
      for (auto spc : spcs->second)
         species.push_back(spc);
    }

    return species;
  }

  map<ResGenus, int> _isospinType;
    int isospinType(ResGenus gen) const {
    auto ptr = _isospinType.find(gen);
    if (ptr == _isospinType.end())
      return (cout << "Isospin not found for particle genus " << gen << endl), 0;
    else
      return ptr->second;
  }

};

}

#endif