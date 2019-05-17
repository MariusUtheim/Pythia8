#ifndef Particle_Widths_H
#define Particle_Widths_H

#include "Pythia8/Info.h"
#include "Pythia8/Interpolator.h"

namespace Pythia8 {

class ParticleWidthEntry {
public:
  Interpolator widths;

  ParticleWidthEntry(Interpolator widthsIn)
    : widths(widthsIn) {}
  ParticleWidthEntry(const ParticleWidthEntry&) = delete;
  ParticleWidthEntry(ParticleWidthEntry&&) = default;

  void addProducts(vector<int> prods, Interpolator brs) {
    branchingRatios.emplace(prods, brs);
  }

  double getWidth(vector<int> prods, double eCM) const {
    auto iter = branchingRatios.find(prods);
    return (iter != branchingRatios.end()) ? iter->second(eCM) * widths(eCM) : 0.;
  }

  double getBR(vector<int> prods, double eCM) const {
    auto iter = branchingRatios.find(prods);
    return (iter != branchingRatios.end()) ? iter->second(eCM) : 0.;
  }

  map<vector<int>, Interpolator> branchingRatios;
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

  double partialWidth(int id, vector<int> prods, double eCM) const;

  double branchingRatio(int id, vector<int> prods, double eCM) const;

  vector<int> getResonances() const;

  vector<vector<int>> getEntry(int id) {
    auto iter = entries.find(id);
    if (iter == entries.end()) return vector<vector<int>>();
    vector<vector<int>> prodss(iter->second.branchingRatios.size());
    for (auto prods : iter->second.branchingRatios)
      prodss.push_back(prods.first);
    return prodss;
  }

  vector<pair<double, vector<int>>> getWeightedProducts(int id, double eCM) const;

  vector<int> pickDecayChannel(int idRes, double eCM);

private:

  Info* infoPtr;

  Rndm* rndmPtr;

  map<int, ParticleWidthEntry> entries;
};

}

#endif