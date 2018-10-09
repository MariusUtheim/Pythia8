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
  
  double dx() const { return (rightSave - leftSave) / ysSave.size(); } 
  double x(int j) const { return leftSave + j * dx(); }

  double operator()(double x) const;

  static const Interpolator Zero;

private:

  double leftSave, rightSave;
  vector<double> ysSave;

};


#endif