#ifndef Mass_Dependent_Width_H
#define Mass_Dependent_Width_H

#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/MassDependentWidth.h"
#include "Pythia8/Interpolator.h"

namespace Pythia8 {

class MassDependentWidth {

public:

  bool readXML(istream& stream);

  double mass(string particle, double eCM) const;

  const Interpolator& getDistribution(string particle) const;

private:

  map<string, Interpolator> massDependentWidths;

};

}

#endif