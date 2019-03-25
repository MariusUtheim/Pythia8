
#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging

int main( int argc, char* argv[] ){
  if (argc != 2) return 1;

  Pythia pythia;
  // Input parameters:
  pythia.readFile(argv[1]);
  // Number of events.
  int nEvent = pythia.mode("Main:numberOfEvents");
  pythia.init();

  // Start generation loop
  for( int iEvent=0; iEvent<nEvent; ++iEvent ){
    cout << "\n\n Begin event no = " << iEvent << endl;

    // Generate next event
    if( !pythia.next() ) {
      if( pythia.info.atEndOfFile() ) break;
      else continue;
    }

    // Debug printout.
    if (iEvent == 25) {
      Vec4 pDiff = -(pythia.process[3].p()+pythia.process[4].p());
      for (int i = 5; i < pythia.process.size(); ++i)
        if (pythia.process[i].isFinal()) pDiff += pythia.process[i].p();
      cout << scientific << setprecision(10)<< " process difference = "
           << pDiff.px() << "  " << pDiff.py() << "  "
           << pDiff.pz() << "  " << pDiff.e() << endl;
      pythia.LHAeventList();
      pythia.process.list();
      pythia.event.list();
    }


  } // end loop over events to generate

  // print cross section, errors
  pythia.stat();

  // Done
  return 0;

}
