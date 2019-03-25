// main program to test new release of Pythia8 in Hidden Valley

#include "Pythia8/Pythia.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace Pythia8;


int main() {


  bool running = false;


  // general settings
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Main:numberOfEvents = 100");
  pythia.readString("Next:numberShowEvent = 5");


  // --- model settings --- //
  pythia.readString("4900001:m0 = 1000");
  pythia.readString("4900001:mWidth = 10");
  pythia.readString("HiddenValley:spinFv = 0");
  pythia.readString("4900001:isResonance = on");
  pythia.readString("4900001:0:bRatio = 1");
  pythia.readString("4900001:0:meMode = 102");

  pythia.readString("4900101:m0 = 10");
  pythia.readString("4900111:m0 = 5");
  pythia.readString("4900113:m0 = 20");
  pythia.readString("4900211:m0 = 5");

  pythia.readString("4900111:tau0 = 150");

  pythia.readString("4900111:0:all on 1.0 102 1 -1");
  pythia.readString("4900113:0:all on 0.999 102 4900111 4900111");
  pythia.readString("4900113:addChannel on 0.001 102 1 -1");
  // ---------------------- //


  // non model-dependent settings
  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("PartonLevel:ISR = on");

  // load processes for emerging jets
  pythia.readString("HiddenValley:gg2DvDvbar = on");
  pythia.readString("HiddenValley:qqbar2DvDvbar = on");

  // run alphaHV
  if (running) {
    pythia.readString("HiddenValley:alphaOrder = 1");
    pythia.readString("HiddenValley:nFlav = 7");
  }
  else pythia.readString("HiddenValley:nFlav = 2");

  // set up HV parton shower
  pythia.readString("HiddenValley:FSR = on");
  pythia.readString("HiddenValley:fragment = on");


  // Initialize
  pythia.init();
  //pythia.particleData.listAll();

  // Begin event loop
  int nEvent = pythia.mode("Main:numberOfEvents");
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
  }

  pythia.stat();

  return 0;

}
