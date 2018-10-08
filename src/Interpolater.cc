
#include "Pythia8/Interpolater.h"

const Interpolater Interpolater::Zero = Interpolater();

Interpolater::Interpolater(istream& stream) {
  double x, y;
  while (stream >> x && stream >> y) {
    xs.push_back(x);
    ys.push_back(y);
  }
}

Interpolater::Interpolater(string path) {
  ifstream stream(path);
  double x, y;
  while (stream >> x && stream >> y) {
    xs.push_back(x);
    ys.push_back(y);
  }
}

double Interpolater::operator()(double x) const
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