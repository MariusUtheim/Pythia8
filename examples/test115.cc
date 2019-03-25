// test115.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Extension of the iTopCopyId and iBotCopyId methods being tested.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    //bool isFlawed = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      for (int j = 0; j < 10; ++j) {
      int iTop1 = pythia.event[i].iTopCopyId(true);
      //int iTop2 = pythia.event[i].iTopCopyId();
      int iBot1 = pythia.event[i].iBotCopyId(true);
      //int iBot2 = pythia.event[i].iBotCopyId();
      }
      /*
      if (iEvent == 0)  cout << " correct copies " << i << "  " << iTop1
        << "  " << iTop2 << "  " << iBot1 << "  " << iBot2 << endl;
      if (iTop2 != iTop1 || iBot2 != iBot1) {
        isFlawed = true;
        cout << " flawed copies " << i << "  " << iTop1 << "  " << iTop2
             << "  " << iBot1 << "  " << iBot2 << endl;
      }
      if (iTop2 > 0 && iTop2 + iBot2 == 0) cout << " error" << endl;
      */
    }
    // if (isFlawed) pythia.event.list();

  // End of event loop. Done.
  }
  return 0;
}
