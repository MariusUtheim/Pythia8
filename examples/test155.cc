#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 1");
  pythia.readString("Beams:frameType = 2");
  pythia.readString("Beams:eA = 4000");
  pythia.readString("Beams:eB = 4000");
  pythia.readString("Bottomonium:O(3S1)[3S1(1)] = 9.28,4.63,3.54");
  pythia.readString("Bottomonium:O(3S1)[3S1(8)] = 0.15,0.045,0.075");
  pythia.readString("Bottomonium:O(3S1)[1S0(8)] = 0.02,0.006,0.01");
  pythia.readString("Bottomonium:O(3S1)[3P0(8)] = 0.02,0.006,0.01");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");
  pythia.readString("Charmonium:gg2ccbar(3S1)[3S1(1)]g = on,off");
  pythia.readString("Charmonium:gg2ccbar(3S1)[3S1(8)]g = on,off");
  pythia.readString("Charmonium:qg2ccbar(3S1)[3S1(8)]q = on,off");
  pythia.readString("Charmonium:qqbar2ccbar(3S1)[3S1(8)]g = on,off");
  pythia.readString("Charmonium:gg2ccbar(3S1)[1S0(8)]g = on,off");
  pythia.readString("Charmonium:qg2ccbar(3S1)[1S0(8)]q = on,off");
  pythia.readString("Charmonium:qqbar2ccbar(3S1)[1S0(8)]g = on,off");
  pythia.readString("Charmonium:gg2ccbar(3S1)[3PJ(8)]g = on,off");
  pythia.readString("Charmonium:qg2ccbar(3S1)[3PJ(8)]q = on,off");
  pythia.readString("Charmonium:qqbar2ccbar(3S1)[3PJ(8)]g = on,off");
  pythia.readString("PDF:useHard = on");
  pythia.readString("PDF:pSet = 13");
  pythia.readString("PartonLevel:all = off");
  pythia.init();
  for (int iEvent = 0; iEvent < 10000; ++iEvent) if (!pythia.next()) continue;
  pythia.stat();
  return 0;
}
