// main16.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates (a) how to collect the analysis code in a separate class
// and (b) how to provide the .cmnd filename on the command line

// Once you have linked the main program you can run it with a command line
// ./main16.exe main16.cmnd > out16

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Put all your own analysis code in the myAnalysis class.

class MyAnalysis {

public:

  // Constructor can be empty.
  //MyAnalysis(): f_eta("strange_eta.dat"),
  //  f_multiplicity("strange_multiplicity.dat") {}
  MyAnalysis() {}

  // Initialization actions.
  void init();

  // Analysis of each new event.
  void analyze(Event& event);

  // Show final results.
  void finish();

private:

  // Declare variables and objects that span init - analyze - finish.
  int  nEvt;
  Hist h_eta, h_multiplicity, h_chgmult;

  //std::ofstream f_eta;
  //std::ofstream f_multiplicity;

  /// below stolen from Rivet

  ///  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
  ///  The Location enum provides a convenient index into the PID.
  enum Location { nj=1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };

  /// Returns everything beyond the 7th digit (e.g. outside the numbering scheme)
  inline int _extraBits(int pid) {
    return abs(pid)/10000000;
  }

  /// Split the PID into constituent integers
  inline unsigned short _digit(Location loc, int pid) {
    //  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj (cf. Location)
    int numerator = (int) std::pow(10.0, (loc-1));
    return (abs(pid)/numerator) % 10;
  }

};

//--------------------------------------------------------------------------

// The initialization code.

void MyAnalysis::init() {

  // Initialize counter for number of events.
  nEvt = 0;

  h_eta.book("strange hadron eta", 30, -6., 6.);
  //h_multiplicity.book("number of strange baryons", 40, 0., 80.);
  h_multiplicity.book("number of strange baryons", 80, -0.5, 79.5);
  h_chgmult.book("number of charged hadrons", 100, -0.5, 799.5);

}

//--------------------------------------------------------------------------

// The event analysis code.

void MyAnalysis::analyze(Event& event) {

  // Increase counter.
  ++nEvt;

  double nb = 0.;
  int nch = 0;

  for(int i=0; i != event.size(); ++i){

    if( !event[i].isFinal() ) continue;
    //    cout<<"is final"<<endl;

    if (event[i].isCharged()) ++nch;

    int pid = event[i].id();

    bool hs = _digit(nq3,pid) == 3 || _digit(nq2,pid) == 3
           || _digit(nq1,pid) == 3;

    if(!hs) continue;

    //cout<<"has strange"<<endl;

    if(_digit(nq3, pid) != 0 ){
    //if(_digit(nq1, pid) != 0 ){
      nb += 1.;
      //cout << " strange baryon id = " << pid << endl;
    }
    h_eta.fill(event[i].eta());
    //f_eta<<event[i].eta()<<endl;
  }

  h_multiplicity.fill(nb);
  //f_multiplicity<<nb<<endl;
  h_chgmult.fill(nch);

}

//--------------------------------------------------------------------------

// The finishing code.

void MyAnalysis::finish() {

  /*
  // Normalize histograms.
  double binFactor = 5. / nEvt;
  yH     *= binFactor;
  etaChg *= binFactor;

  // Print histograms.
  cout << brH << yH << etaChg << mult;
*/

  cout << h_eta << endl<<h_multiplicity<<h_chgmult<<endl ;

}

//==========================================================================

// You should not need to touch the main program: its actions are
// determined by the .cmnd file and the rest belongs in MyAnalysis.

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file name and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided file name corresponds to an existing file.
  ifstream is(argv[1]);
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Confirm that external file will be used for settings..
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;

  // Declare generator. Read in commands from external file.
  Pythia pythia;
  pythia.readFile(argv[1]);

  // Initialization.
  pythia.init();

  // Declare user analysis class. Do initialization part of it.
  MyAnalysis myAnalysis;
  myAnalysis.init();

  // Read in number of event and maximal number of aborts.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");
  bool hasPL = pythia.flag("PartonLevel:all");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // User Analysis of current event.
    myAnalysis.analyze( (hasPL ? pythia.event : pythia.process) );

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();

  // User finishing.
  myAnalysis.finish();

  // Done.
  return 0;
}
