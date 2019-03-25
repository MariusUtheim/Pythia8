// main45.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Stephen Mrenna, mrenna@fnal.gov
// This program illustrates how to call Rivet directly
// from Pythia8 without introducing fifos or hepmc files.

#include "Pythia8/Pythia.h"
#include <iostream>

using namespace Pythia8;


int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide one input (parameter) file name. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided input name corresponds to an existing file.
  ifstream is(argv[1]);
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Generator.
  Pythia pythia;


  //  Event &event = pythia.event;
  // Read in commands from external file.
  pythia.readFile(argv[1]);


  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");

  // Initialization.
  pythia.init();

  Hist _h1000993("/RHadron/93 Mass - offset [GeV]",100,1.0,3.0);
  Hist _h1009113("/RHadron/113 Mass - offset [GeV]",100,1.0,3.0);
  Hist _h1009213("/RHadron/213 Mass - offset [GeV]",100,1.0,3.0);
  Hist _h1009223("/RHadron/223 Mass - offset [GeV]",100,1.0,3.0);
  Hist _h1009323("/RHadron/323 Mass - offset [GeV]",100,1.0,3.0);
  Hist _h1092214("/RHadron/2214 Mass - offset  [GeV]",100,1.0,3.0);
  Hist _h1092224("/RHadron/2224 Mass - offset  [GeV]",100,1.0,3.0);

  double mOffset = pythia.particleData.m0(1000021);
  cout << " 1000993  " << pythia.particleData.m0(1000993) << endl;
  cout << " 1009113  " << pythia.particleData.m0(1009113) << endl;
  cout << " 1009213  " << pythia.particleData.m0(1009213) << endl;
  cout << " 1009223  " << pythia.particleData.m0(1009223) << endl;
  cout << " 1009323  " << pythia.particleData.m0(1009323) << endl;
  cout << " 1092214  " << pythia.particleData.m0(1092214) << endl;
  cout << " 1092224  " << pythia.particleData.m0(1092224) << endl;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate event.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
        cout << " Aborted since reached end of Les Houches Event File\n";
        break;
      }

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
    for (int i = 0; i < pythia.event.size(); ++i) {
      int idLocal = pythia.event[i].idAbs();
      double mFill = pythia.event[i].m() - mOffset;
      switch (idLocal) {
      case 1000993:
	_h1000993.fill( mFill );
	break;
      case 1009113:
	_h1009113.fill( mFill );
	break;
      case 1009213:
	_h1009213.fill( mFill );
	break;
      case 1009223:
	_h1009223.fill( mFill );
	break;
      case 1009323:
	_h1009323.fill( mFill );
	break;
      case 1092214:
	_h1092214.fill( mFill );
	break;
      case 1092224:
	_h1092224.fill( mFill );
	break;
      default:
	break;
      }
    }

  // End of event loop. Statistics.
  }
  pythia.stat();
//TEST
  cout << "finalizing" << endl;
  string rivetOutput=pythia.word("Main:spareWord2");

  std::cout << _h1000993;
  std::cout << _h1009113;
  std::cout << _h1009213;
  std::cout << _h1009223;
  std::cout << _h1009323;
  std::cout << _h1092214;
  std::cout << _h1092224;


  // Done.
  return 0;
}
