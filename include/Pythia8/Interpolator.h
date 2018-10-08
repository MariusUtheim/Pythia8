#ifndef Pythia8_Interpolator_H
#define Pythia8_Interpolator_H

#include "PythiaStdlib.h"
#include <istream>

using namespace Pythia8;

class Interpolator {
public:

  Interpolator(string path);
  Interpolator(istream& stream);

  double operator()(double x) const;

  static const Interpolator Zero;

private:
  
  Interpolator() { xs = { 0. }; ys = { 0. }; }

  vector<double> xs, ys;

};


#endif