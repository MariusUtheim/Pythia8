#ifndef Pythia8_Interpolator_H
#define Pythia8_Interpolator_H

#include "PythiaStdlib.h"
#include <istream>

using namespace Pythia8;

class Interpolator {
public:

  Interpolator(double leftIn, double rightIn, vector<double> ysIn);

  double left() { return leftSave; }
  double right() { return rightSave; }
  
  double dx() const { return (rightSave - leftSave) / (ysSave.size() - 1); } 
  double x(int j) const { return leftSave + j * dx(); }
  const vector<double>& data() const { return ysSave; }

  double operator()(double x) const;

  static const Interpolator Zero;

private:

  const double leftSave, rightSave;
  const vector<double> ysSave;

};


#endif