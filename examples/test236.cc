#include "Pythia8/Pythia.h"
using namespace Pythia8;

// Test the MPI mass changes.
void mpi() {
  // Configure Pythia, artificially inflate quarkonia.
  Pythia pythia("", false);
  pythia.readString("Beams:eCM = 8000.");
  //pythia.readString("Print:quiet = on");
  pythia.readString("SoftQCD:all = on");
  pythia.readString("Charmonium:O(3S1)[3S1(1)] = 1e5,1e5");
  pythia.init();

  // Check isOnium returns correct result.
  ParticleData &pd = pythia.particleData;
  int id(1);
  while (id) {
    if (pd.isOnium(id))
      cout << setw(15) << pd.name(id) << setw(10) << id << "\n";
    id = pd.nextId(id);
  }

  // Begin event loop. Generate event. Skip if error. List first one.
  Event &event = pythia.process;
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    for (int iPrt = 0; iPrt < event.size(); ++iPrt) {
      if (pd.isOnium(event[iPrt].id()))
	cout << setw(10) << event[iPrt].id() << scientific
	     << setprecision(5) << setw(15) << event[iPrt].mCalc() << "\n";
    }
  }
}


// Test the hard process changes.
void hard() {
  // Configure Pythia, artificially inflate quarkonia.
  Pythia pythia("", false);
  pythia.readString("Beams:eCM = 8000.");
  //pythia.readString("Print:quiet = on");
  pythia.readString("Charmonium:gg2ccbar(3S1)[3S1(1)]g = on, off");
  //pythia.readString("PhaseSpace:minWidthBreitWigners = 1e-6");
  pythia.init();

  // Begin event loop. Generate event. Skip if error. List first one.
  ParticleData &pd = pythia.particleData;
  Event &event = pythia.process;
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    for (int iPrt = 0; iPrt < event.size(); ++iPrt) {
      if (pd.isOnium(event[iPrt].id()))
	cout << setw(10) << event[iPrt].id() << scientific
	     << setprecision(5) << setw(15) << event[iPrt].mCalc() << "\n";
    }
  }
}

// The main program.
int main() {

  // Run the tests.
  mpi();
  hard();

  // Finalize.
  return 0;
}
