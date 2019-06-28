// Interpolator.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Interpolator class.

#ifndef Pythia8_Interpolator_H
#define Pythia8_Interpolator_H

#include "PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================
// @TBD Maybe move this to Basics.h?
// Interpolator class.
// Used to interpolate between values in linearly spaced data.

class Interpolator {
public:

  // Constructor.
  Interpolator(double leftIn, double rightIn, vector<double> ysIn)
    : leftSave(leftIn), rightSave(rightIn), ysSave(ysIn) { }

  // Function to get y-values of interpolation data. 
  const vector<double>& data() const { return ysSave; }

  // x-values are linearly spaced on the interpolation region
  double left()  const { return leftSave; }
  double right() const { return rightSave; }

  // Operator to get interpolated value at the specified point
  double operator()(double x) const;

private:

  // Data members
  double leftSave, rightSave;
  vector<double> ysSave;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Interpolator_H
