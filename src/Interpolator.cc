
#include "Pythia8/Interpolator.h"

const Interpolator Interpolator::Zero(0., 1., { 0., 0. });

Interpolator::Interpolator(double leftIn, double rightIn, vector<double> ysIn)
  : leftSave(leftIn), rightSave(rightIn), ysSave(ysIn) {
  if (ysIn.size() <= 1)
    cout << "WARNING: Interpolator has too few entries! " << ysIn.size() << endl;
}

double Interpolator::operator()(double xIn) const {

  double t = (xIn - leftSave) / (rightSave - leftSave);
  int j = (int)floor(t * (ysSave.size() - 1));

  if (j < 0)
    return ysSave[0];
  else if (j >= ysSave.size() - 1)
    return ysSave[ysSave.size() - 1];
  else {
    double s = (xIn - this->x(j)) / dx();
    return (1 - s) * ysSave[j] + s * ysSave[j + 1];
  }
}