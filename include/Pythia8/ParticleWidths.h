#ifndef Particle_Widths_H
#define Particle_Widths_H

#include "Pythia8/Info.h"
#include "Pythia8/Interpolator.h"

namespace Pythia8 {

class ParticleWidthEntry {
public:

  ParticleWidthEntry(Interpolator widthsIn)
    : widths(widthsIn) {}
  ParticleWidthEntry(const ParticleWidthEntry&) = delete;
  ParticleWidthEntry(ParticleWidthEntry&&) = default;

  void addProducts(pair<int, int> prods, Interpolator brs) {
    branchingRatios.emplace(prods, brs);
  }

  double getWidth(pair<int, int> prods, double eCM) const {
    auto iter = branchingRatios.find(prods);
    return (iter != branchingRatios.end()) ? iter->second(eCM) * widths(eCM) : 0.;
  }

  double getBR(pair<int, int> prods, double eCM) const {
    auto iter = branchingRatios.find(prods);
    return (iter != branchingRatios.end()) ? iter->second(eCM) : 0.;
  }

  double m0;
  Interpolator widths;
  map<pair<int, int>, Interpolator> branchingRatios;
};

//TS?? HadronWidths better name?
//Then also ParticleWidths.xml -> HadronWidthData.xml
class ParticleWidths {

public:

  bool init(Info* infoPtrIn, Rndm* rndmPtrIn, string path) 
  {
    infoPtr = infoPtrIn;
    rndmPtr = rndmPtrIn;
    ifstream stream(path);
    if (!stream.is_open()) {
      infoPtr->errorMsg( "Warning in ParticleWidths::init: "
          "unable to open file");
      return false;
    }
    return readXML(stream);
  }

  bool readXML(istream& stream);

  double width(int id, double eCM) const;

  double partialWidth(int id, pair<int, int> prods, double eCM) const;

  double branchingRatio(int id, pair<int, int> prods, double eCM) const;

  vector<int> getResonances() const;

  // @TODO Make this work with antiparticles
  vector<pair<int, int>> getEntry(int id) {
    auto iter = entries.find(id);
    if (iter == entries.end()) return vector<pair<int, int>>();
    vector<pair<int, int>> prodss(iter->second.branchingRatios.size());
    for (auto prods : iter->second.branchingRatios)
      prodss.push_back(prods.first);
    return prodss;
  }

  // @TODO: Make this work with antiparticles
  bool hasData(int id) const {
    auto iter = entries.find(id);
    return iter != entries.end();
  }

  vector<pair<double, pair<int, int>>> getWeightedProducts(int id, double eCM) const;

  pair<int, int> pickDecayChannel(int idRes, double eCM);

  double pickMass(int idRes, double eCM, double mB, int lType = 1);

  pair<double, double> pickMass2(int id1, int id2, double eCM, int lType = 1);

  bool pickDecay(int idDec, double m, int& id1Out, int& id2Out,
    double& m1Out, double& m2Out);

private:

  Info* infoPtr;

  Rndm* rndmPtr;

  map<int, ParticleWidthEntry> entries;
};

}

#endif