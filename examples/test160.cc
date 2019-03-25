// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString(" HardQCD:all = on");
  pythia.readString("    PhaseSpace:pTHatMin = 20.");
  pythia.readString("  PartonLevel:MPI = off");

  // Play with words and vectors of words.
  pythia.readString("PDF:piSet = {string with blanks}");
  pythia.readString("Bottomonium:O(3S1)[3S1(1)] = { 5.52 , 27.7, 2.3 }");
  pythia.settings.addWVec("sentence", vector<string>(1," "));
  pythia.readString("sentence = { A sentence, with=a=beginning, and an end. }");
  pythia.settings.addWVec("  UncertaintyVariation", vector<string>(1," "));
  pythia.readString("  UncertaintyVariation = { A fsr:muR 0.5 isr:muR 0.5 ,");
  pythia.readString(" C fsr:cNS 3.0}");
  pythia.settings.addWVec("input2", vector<string>(1," "));
  pythia.readFile("test160.cmnd");

  // Initialize for printout.
  pythia.init();

  // Recover input.
  vector<string> sentence = pythia.settings.wvec("sentence");
  for (int i = 0; i < int(sentence.size()); ++i)
    cout << " part " << i << " is \"" << sentence[i] << "\"" << endl;
  vector<string> uncert = pythia.settings.wvec("UncertaintyVariation");
  for (int i = 0; i < int(uncert.size()); ++i)
    cout << " part " << i << " is \"" << uncert[i] << "\"" << endl;

  return 0;
}
