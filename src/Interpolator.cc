// Interpolator.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Implementation of the interpolation function used by th eInterpolator class.

#include "Pythia8/Interpolator.h"

namespace Pythia8 {

//==========================================================================

// Interpolator class.
// Used to interpolate between values in linearly spaced data.

//--------------------------------------------------------------------------

// Operator to get interpolated value at the specified point

double Interpolator::operator()(double xIn) const {
  // @TODO Decide what happens if out of range
  double t = (xIn - leftSave) / (rightSave - leftSave);
  int lastIdx = ysSave.size() - 1;
  int j = (int)floor(t * lastIdx);
  double dx = (rightSave - leftSave) / (ysSave.size() - 1);

  if (j < 0)
    return ysSave[0];
  else if (j >= lastIdx)
    return ysSave[lastIdx];
  else {
    double s = (xIn - (leftSave + j * dx)) / dx;
    return (1 - s) * ysSave[j] + s * ysSave[j + 1];
  }
}

//==========================================================================

} // end namespace Pythia8
