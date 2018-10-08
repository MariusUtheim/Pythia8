
#include "Pythia8/Interpolator.h"

const Interpolator Interpolator::Zero = Interpolator();

Interpolator::Interpolator(istream& stream) {
  double x, y;
  while (stream >> x && stream >> y) {
    xs.push_back(x);
    ys.push_back(y);
  }
}

Interpolator::Interpolator(string path) {
  ifstream stream(path);
  double x, y;
  while (stream >> x && stream >> y) {
    xs.push_back(x);
    ys.push_back(y);
  }
}

double Interpolator::operator()(double x) const
{
  if (x < xs[0])
    return ys[0];

  for (int i = 0; i < xs.size() - 1; ++i) {
    if (x >= xs[i] && x < xs[i + 1]) {
      double t = (x - xs[i]) / (xs[i + 1] - xs[i]);
      return (1 - t) * ys[i] + t * ys[i + 1];
    }
  }

  return ys[xs.size() - 1];
}