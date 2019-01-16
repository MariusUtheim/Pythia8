#ifndef Mass_Dependent_Width_H
#define Mass_Dependent_Width_H

#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Interpolator.h"

namespace Pythia8 {

class MassDependentWidth {

public:

  bool init(string path) {
    ifstream stream(path);
    if (!stream.is_open()) return false;
    return readXML(stream);
  }

  bool readXML(istream& stream);

  double width(string particle, double eCM) const;

  double branchingRatio(string particle, string products, double eCM) const;

  const Interpolator& getDistribution(string particle) const;

  const Interpolator& getBranchingRatios(string particle, string brs) const;

  const vector<string>& getProducts(string particle) const;

private:

  map<string, vector<string>> productList;

  map<string, Interpolator> massDependentWidths;

  map<pair<string, string>, Interpolator> branchingRatios;
};

}

#endif