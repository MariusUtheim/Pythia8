
#include "Pythia8/Interpolator.h"

const Interpolator Interpolator::Zero = Interpolator(0., 1., { 0. });

Interpolator::Interpolator(double leftIn, double rightIn, vector<double> ysIn)
  : leftSave(leftIn), rightSave(rightIn), ysSave(ysIn) { }

double Interpolator::operator()(double x) const {

  double t = (x - leftSave) / (rightSave - leftSave);
  if (t < 0)
    return ysSave[0];
  else if (t >= 1)
    return ysSave[ysSave.size() - 1];
  else {
    int j = (int)(t * ysSave.size() - 1);
    double s = (x - this->x(j)) / dx();
    return (1 - s) * ysSave[j] + s * ysSave[j + 1];
  }
}