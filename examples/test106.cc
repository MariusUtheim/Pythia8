// test106.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Check that tunes are set correctly.
#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Loop over tunes.
  for (int ppTune = 1; ppTune < 33; ++ppTune) {

    // Generator. Empty process. Tune. Initialize.
    Pythia pythia;
    pythia.readString("ProcessLevel:all = off");
    pythia.settings.mode("Tune:pp", ppTune);
    pythia.init();

  // Done.
  }
  return 0;
}
