#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {
  Pythia pythia;
  pythia.readString("Next:showScaleAndVertex = on");
  pythia.readString("Next:numberShowEvent = 1");

  pythia.readString("HiddenValley:Ngauge = 1");
  pythia.readString("HiddenValley:ffbar2Zv = on");
  pythia.readString("HiddenValley:FSR = on");
  pythia.readString("HiddenValley:alphaFSR = 0.2");

  pythia.readString("4900101:m0 = 0.4");
  pythia.readString("4900023:m0 = 1000.");
  pythia.readString("4900023:mMin = 500.");
  pythia.readString("4900023:mWidth = 31.3");
  pythia.readString("4900023:onMode = off");
  pythia.readString("4900023:OnIfMatch = 4900101 -4900101");

  pythia.readString("4900022:m0 = 0.1");
  pythia.readString("4900022:tau0 = 1.");
  pythia.readString("4900022:oneChannel = 1 1 91 11 -11");

  pythia.init();
  double Ngenerate = 1;
  for(int iEvent = 1; iEvent <= Ngenerate; ++iEvent){
	if( !pythia.next() ) continue;
	if( iEvent > Ngenerate ) break;
  }//loop over events

  pythia.stat();
  return 0;
}
