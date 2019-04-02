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

  void addProducts(vector<int> prods, Interpolator brs) {
    branchingRatios.emplace(prods, brs);
  }

  double getBR(vector<int> prods, double eCM) const {
    auto entry = branchingRatios.find(prods);
    if (entry == branchingRatios.end())
      return 0;
    else
      return entry->second(eCM);
  }

  map<vector<int>, Interpolator> branchingRatios;
};

class ParticleWidths {

public:

  bool init(Info* infoPtrIn, string path) 
  {
    infoPtr = infoPtrIn;
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

  double branchingRatio(int id, vector<int> prods, double eCM) const;

  vector<int> getResonances() const;

  vector<pair<double, vector<int>>> getWeightedProducts(int id, double eCM) const;

private:

  Info* infoPtr;

  map<int, ParticleWidthEntry> entries;
};

}

#endif