// test119.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how to write a runtime Les Houches interface.
// Generates a fixed configuration over and over again.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Write own derived LHAup class.

class LHAupDY : public LHAup {

public:

  // Constructor and destructor trivial.
  LHAupDY() {}
  ~LHAupDY() {}

  // Initialization of incoming beams: pp collisions at 8 TeV.
  // Also specify process, here with fictitious cross section = 1.
  virtual bool setInit() {
    eCM = 8000.;
    setBeamA( 2212, 0.5 * eCM);
    setBeamB( 2212, 0.5 * eCM);
    addProcess( 1, 1., 0., 1.);
    return true;
  }

  // Set the incoming event. Here same event over and over again.
  virtual bool setEvent(int = 0) {

    // Specify incoming flavour and x fractions.
    int idIn  = 2;
    double x1 = 0.05;
    double x2 = 0.003;

    // Specify hard process, including its ficitition alpha_em and alpha_s.
    double scale   = sqrt(x1 * x2) * eCM;
    double alphaEM = 1./128.;
    double alphaS  = 0.12;
    setProcess( 1, 1., scale, alphaEM, alphaS);

    // Store incoming and produced particles.
    addParticle(  idIn, -1, 0, 0, 101, 0,  0., 0.,  0.5 * x1 * eCM,
       0.5 * x1 * eCM, 0., 0., 9., scale);
    addParticle( -idIn, -1, 0, 0, 0, 101,  0., 0., -0.5 * x2 * eCM,
       0.5 * x2 * eCM, 0., 0., 9., scale);
    addParticle( 23, 1, 1, 2, 0, 0,  0., 0., 0.5 * (x1 - x2) * eCM,
       0.5 * (x1 + x2) * eCM, scale, 0., 9., scale);

    // Set info on incoming beams (not required) and done.
    setPdf( idIn, -idIn, x1, x2, scale, 1., 1., true);
    return true;
  }

  // Stored quantities.
  double eCM;

};

//--------------------------------------------------------------------------

int main() {

  // Generator.
  Pythia pythia;

  // Create instance of runtime Les Houches generator.
  LHAup* lhaup = new LHAupDY();
  pythia.setLHAupPtr( lhaup);

  // Switch off Z0 decay products.
  pythia.readString("23:mayDecay = off");

  // Initialize.
  pythia.readString("Beams:frameType = 5");
  pythia.init();

  // Book histogram.
  Hist nCharged( "charged multiplicity", 100, -1., 399.);

  // Number of events to generate. Allow for a few faulty events.
  int nEvent = 100;
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events, and check whether generation failed.
    // First few failures written off as "acceptable" errors, then quit.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      break;
    }

    // Sum up final charged multiplicity and fill in histogram.
    int nChg = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
      ++nChg;
    nCharged.fill(nChg);

  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();
  cout << nCharged;

  // Done.
  return 0;
}
