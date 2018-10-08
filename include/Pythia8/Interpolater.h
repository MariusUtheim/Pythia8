#ifndef Pythia8_Interpolater_H
#define Pythia8_Interpolater_H

#include "PythiaStdlib.h"
#include <istream>

using namespace Pythia8;

class Interpolater {
public:

  Interpolater(string path);
  Interpolater(istream& stream);

  double operator()(double x) const;

  static const Interpolater Zero;

private:
  
  Interpolater() { xs = { 0. }; ys = { 0. }; }

  vector<double> xs, ys;

};


#endif