// Debug program from Ryosuke Sato.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {
  int nEvent = 100000000;
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("ProcessLevel:all = off");
  //pythia.readString("ParticleDecays:FSRinDecays = off");
  pythia.init();

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    event.reset();
    // Upsilon(1S) at rest
    event.append( 553, 1, 0, 0, 0., 0., 0., 9.46030, 9.46030);
    if (!pythia.next()) --iEvent;
  }

  return 0;
}
