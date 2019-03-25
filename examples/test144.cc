// test144.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Debug high cross section for Andy Buckley.
// Reason related to long, dominant low-mass tail badly caught.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Key settings to be used in the main program.
  int nEvent   = 1000;
  int nAbort   = 3;
  double eCM   = 8000.;

  // Generator. Shorthand for the event.
  Pythia pythia;

  // Set up beams: p p is default so only need set energy.
  pythia.settings.parm("Beams:eCM", eCM);

  // SUSY production limited to g g -> ~c_L ~cbar_L.
  pythia.readString("SUSY:gg2squarkantisquark = on");
  pythia.readString("SUSY:idA =  1000004");
  pythia.readString("SUSY:idB = -1000004");

  // Troublemaker SLHA file.
  pythia.readString("SLHA:file = test144.spc");

  // Mass range.
  pythia.readString("1000004:mMin = 200.");
  //pythia.readString("1000004:mMax = 3100.");

  // Various checks.
  pythia.readString("PhaseSpace:increaseMaximum = on");
  pythia.readString("PhaseSpace:showViolation = on");

  // Switch off key components.
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Init:showOneParticleData = 1000004");

  // Initialize.
  pythia.init();

  // Histograms.
  Hist mHat( "mHat", 100, 0., 10000.);
  Hist m3( "m3", 100, 0., 5000.);
  Hist m4( "m4", 100, 0., 5000.);
  Hist mDif( "mDif", 100, 0., 1000.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    //cout << " Begin event no = " << iEvent << endl;

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Study event.
    mHat.fill( pythia.info.mHat() );
    m3.fill( pythia.process[5].m() );
    m4.fill( pythia.process[6].m() );
    mDif.fill( pythia.info.mHat() - pythia.process[5].m()
       - pythia.process[6].m() );

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();
  cout << mHat << m3 << m4 << mDif;

  // Done.
  return 0;
}
