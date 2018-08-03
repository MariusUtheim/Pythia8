// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include <fstream>
using namespace Pythia8;


struct SimpleEventListing {
  int no;
  int mother;
  bool isFinal;
  int id;
  double e, px, py, pz;

  SimpleEventListing(int no, int mother, int isFinal, int id, double e, double px, double py, double pz) :
    no(no), mother(mother), isFinal(isFinal), id(id), e(e), px(px), py(py), pz(pz) {}

  void write(char* buffer) {
    char *c = reinterpret_cast<char *>(this);
    for (int i = 0; i < sizeof(SimpleEventListing); i++)
      buffer[i] = c[i];
  }
};



int main() {
  Pythia pythia;
  pythia.readFile("mymain01.cmnd");
  pythia.init();
  
  while (!pythia.next());

  Event ev = pythia.event[0];

  std::ofstream f("pythia.bin", std::ios::binary);
  for (int i = 0; i < ev.size(); i++)
  {
    Particle p = ev[i];
    SimpleEventListing s(i, p.mother1(), p.isFinal(), p.id(), p.e(), p.px(), p.py(), p.pz());
    f.write(&s, sizeof(SimpleEventListing));
  }

  f.close();
  
  pythia.stat();
  return 0;
}
